{-# LANGUAGE ScopedTypeVariables #-}
module Jack2
  (jack, zonal, schur)
  where
import           Control.Lens                             hiding (iconcatMap)
import           Control.Monad                            (when)
import           Data.Array.IO
import           Data.List.Index                          (iconcatMap)
import           Data.Maybe
import           Math.Combinat.Partitions.Integer.IntList (_dualPartition,
                                                           _isPartition)
import           Numeric.SpecFunctions                    (factorial)

_ij :: [Int] -> ([Int],[Int])
_ij lambda =
  (
    iconcatMap (\i a ->  replicate a (i+1)) lambda,
    concatMap (\a -> [1 .. a]) (filter (>0) lambda)
  )

_convParts :: Num b => [Int] -> ([b],[b])
_convParts lambda =
  (map fromIntegral lambda, map fromIntegral (_dualPartition lambda))

hookLengths :: Fractional a => [Int] -> a -> [a]
hookLengths lambda alpha = upper ++ lower
  where
  (i,j) = _ij lambda
  (lambda', lambdaConj') = _convParts lambda
  upper = zipWith (fup lambdaConj' lambda') i j
    where
    fup x y ii jj =
      x!!(jj-1) - fromIntegral ii + alpha*(y!!(ii-1) - fromIntegral jj + 1)
  lower = zipWith (flow lambdaConj' lambda') i j
    where
    flow x y ii jj =
      x!!(jj-1) - fromIntegral ii + 1 + alpha*(y!!(ii-1) - fromIntegral jj)

_B :: Fractional a => [Int] -> [Int] -> [Int] -> a -> a
_B nu lambda mu alpha = if all (== 0) nu then 1 else product hookLngths
  where
  (i,j) = _ij nu
  (nu', nuConj) = _convParts nu
  lambdaConj = _dualPartition lambda ++ repeat 0
  muConj = _dualPartition mu ++ repeat 0
  hookLngths = zipWith (f lambdaConj muConj) i j
    where
    f a b ii jj = if a!!(jj-1) == b!!(jj-1)
      then
        fup nuConj nu' ii jj
      else
        flow nuConj nu' ii jj
      where
      fup x y k l =
        x!!(l-1) - fromIntegral k + alpha*(y!!(k-1) - fromIntegral l + 1)
      flow x y k l =
        x!!(l-1) - fromIntegral k + 1 + alpha*(y!!(k-1) - fromIntegral l)

_beta :: Fractional a => [Int] -> [Int] -> a -> a
_beta lambda mu alpha = _B lambda lambda mu alpha / _B mu lambda mu alpha

_N :: [Int] -> [Int] -> Int
_N lambda mu = sum $ zipWith (*) mu prods
  where
  prods = map (\i -> product $ drop i (map (+1) lambda)) [1 .. length lambda]

jack :: forall a. (Fractional a, Ord a) => [a] -> [Int] -> a -> IO a
jack x lambda alpha =
  case _isPartition lambda && alpha > 0 of
    False -> if _isPartition lambda
      then error "alpha must be strictly positive"
      else error "lambda is not a valid integer partition"
    True -> do
      arr0 <- newArray ((1,1), (_N lambda lambda, length x)) Nothing
      jac (length x) 0 lambda lambda arr0 Nothing
        where
        theproduct :: Int -> a
        theproduct nu0 = if nu0 <= 1
          then 1
          else product $ map (\i -> alpha * fromIntegral i + 1) [1 .. nu0-1]
        jac :: Int -> Int -> [Int] -> [Int] -> IOArray (Int,Int) (Maybe a)
               -> Maybe a -> IO a
        jac m k mu nu arr elt
          | nu!!0 == 0 || m == 0 = return 1
          | length nu > m && nu!!m > 0 = return 0
          | m == 1 = return $ x!!0^(nu!!0) * theproduct (nu!!0)
          | k == 0 && isJust elt = do
            e <- readArray arr (_N lambda nu, m)
            return $ fromJust e
          | otherwise = do
            e <- readArray arr (_N lambda nu, m-1)
            jck <- jac (m-1) 0 nu nu arr e
            let s = jck * _beta mu nu alpha * x!!(m-1)^(sum mu - sum nu)
            ss <- go s (max 1 k)
            when (k == 0) $ writeArray arr (_N lambda nu, m) (Just ss)
            return ss
            where
            go :: a -> Int -> IO a
            go ss ii
              | length nu < ii || nu!!(ii-1) == 0 = return ss
              | otherwise =
                if length nu == ii && nu!!(ii-1)>0 || nu!!(ii-1)>nu!!ii
                  then do
                    let nu' = (element (ii-1) .~ nu!!(ii-1)-1) nu
                    if nu!!(ii-1) > 1
                      then do
                        e <- readArray arr (_N lambda nu', m)
                        jck <- jac m ii mu nu' arr e
                        go (ss+jck) (ii+1)
                      else
                        if nu'!!0 == 0
                          then do
                            let jck' = _beta mu nu' alpha * x!!(m-1)^ sum mu
                            go (ss+jck') (ii+1)
                          else do
                            e <- readArray arr (_N lambda nu', m-1)
                            jck <- jac (m-1) 0 nu' nu' arr e
                            let jck' = jck * _beta mu nu' alpha *
                                         x!!(m-1)^(sum mu - sum nu')
                            go (ss+jck') (ii+1)
                  else
                    go ss (ii+1)

zonal :: (Fractional a, Ord a) => [a] -> [Int] -> IO a
zonal x lambda = do
  let k = sum lambda
      jlambda = product (hookLengths lambda 2)
      c = 2^k * realToFrac (factorial k) / jlambda
  jck <- jack x lambda 2
  return $ c * jck

schur :: forall a. Fractional a => [a] -> [Int] -> IO a
schur x lambda =
  case _isPartition lambda of
    False -> error "lambda is not a valid integer partition"
    True -> do
      arr0 <- newArray ((1,1), (_N lambda lambda, length x)) Nothing
      sch (length x) 1 lambda arr0 Nothing
        where
        sch :: Int -> Int -> [Int] -> IOArray (Int,Int) (Maybe a)
               -> Maybe a -> IO a
        sch m k nu arr elt
          | nu!!0 == 0 || m == 0 = return 1
          | length nu > m && nu!!m > 0 = return 0
          | m == 1 = return $ x!!0^(nu!!0)
          | isJust elt = do
            e <- readArray arr (_N lambda nu, m)
            return $ fromJust e
          | otherwise = do
            e <- readArray arr (_N lambda nu, m-1)
            s <- sch (m-1) 1 nu arr e
            ss <- go s k
            when (k == 1) $ writeArray arr (_N lambda nu, m) (Just ss)
            return ss
            where
            go :: Fractional a => a -> Int -> IO a
            go ss ii
              | length nu < ii || nu!!(ii-1) == 0 = return ss
              | otherwise =
                if length nu == ii && nu!!(ii-1) > 0 || nu!!(ii-1) > nu!!ii
                  then do
                    let nu' = (element (ii-1) .~ nu!!(ii-1)-1) nu
                    if nu!!(ii-1) > 1
                      then do
                        e <- readArray arr (_N lambda nu', m)
                        schr <- sch m ii nu' arr e
                        go (ss + x!!(m-1)*schr) (ii+1)
                      else
                        if nu'!!0 == 0
                          then
                             go (ss + x!!(m-1)) (ii+1)
                          else do
                            e <- readArray arr (_N lambda nu', m-1)
                            schr <- sch (m-1) 1 nu' arr e
                            go (ss + x!!(m-1)*schr) (ii+1)
                  else
                    go ss (ii+1)
