module GPochhammer where
import           Jack4                                    (zonal)
import           Math.Combinat.Partitions.Integer.IntList (_allPartitions)

gpochhammer :: Fractional a => a -> [Int] -> a -> a
gpochhammer a kappa alpha =
  product $
    map (\i -> product $
                 map (\j -> a - (fromIntegral i - 1)/alpha + fromIntegral j -1)
                     [1 .. kappa !! (i-1)])
        [1 .. length kappa]

hcoeff :: Fractional a => [a] -> [a] -> [Int] -> a -> a
hcoeff a b kappa alpha =
  numerator / denominator / fromIntegral (factorial (sum kappa))
  where
    factorial n = product [1 .. n]
    numerator = product $ map (\x -> gpochhammer x kappa alpha) a
    denominator = product $ map (\x -> gpochhammer x kappa alpha) b

testHypergeo :: Double
testHypergeo =
  let a = [2,3] in
  let b = [4] in
  let coeff kappa = hcoeff a b kappa 2 in
  let kappas = [[], [1], [1,1], [2]] in
  let x = [5,6] in
  sum $ map (\kappa -> coeff kappa * zonal x kappa) kappas

hypergeoPQ :: (Fractional a, Ord a) => Int -> [a] -> [a] -> [a] -> a
hypergeoPQ m a b x =
  sum $ map (\kappa -> coeff kappa * zonal x kappa) kappas
  where
  kappas = filter (\kap -> length kap <= length x) (_allPartitions m)
  coeff kappa = hcoeff a b kappa 2
