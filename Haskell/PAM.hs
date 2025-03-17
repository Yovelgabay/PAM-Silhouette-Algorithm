import Data.List
import Data.List.Split
import Data.Ord
import Data.Time.Clock
import Numeric.LinearAlgebra
import System.IO
import Prelude hiding ((<>))

-- | Read a matrix from a text file and return it as a list of lists.
readDistanceMatrix :: (Floating a, Read a) => FilePath -> String -> IO [[a]]
readDistanceMatrix path delim = do
  file <- openFile path ReadMode
  text <- hGetContents file
  return $ map (map read . splitOn delim) (lines text)

-- | Extract a submatrix from a matrix based on a list of indices.
extractSubmatrix :: Matrix Double -> [Int] -> Matrix Double
extractSubmatrix matrix indices = (matrix ? indices) ¿ indices

-- | Calculate the sum of distances for a given list of points from a specific point.
sumDistances :: Matrix Double -> Int -> [Int] -> Double
sumDistances distanceMatrix p points =
  sumElements (fromList [distanceMatrix `atIndex` (p, i) | i <- points])

-- | Find the point with the minimum sum of distances in a list of points.
findNewMedoid :: Matrix Double -> [Int] -> Int
findNewMedoid distanceMatrix points =
  points !! minIndex (fromList [sumDistances distanceMatrix p points | p <- points])

-- | Get all points assigned to a specific medoid.
getClusterPoints :: [Int] -> Int -> [Int]
getClusterPoints clusters medoidVal = elemIndices medoidVal clusters

-- | Find the optimal medoid for each cluster.
findMedoid :: Int -> Matrix Double -> [Int] -> [Int] -> [Int]
findMedoid n distanceMatrix clusters medoids = findMedoids 0 medoids
  where
    findMedoids m updatedMedoids
      | m >= length medoids = updatedMedoids
      | otherwise = findMedoids (m + 1) (replaceAtIndex m newMedoid updatedMedoids)
      where
        points = getClusterPoints clusters (updatedMedoids !! m)
        newMedoid = findNewMedoid distanceMatrix points
    replaceAtIndex i newVal updatedM = take i updatedM ++ [newVal] ++ drop (i + 1) updatedM

-- | Assign each point to the nearest medoid.
assignPointsToClusters :: Matrix Double -> [Int] -> Int -> [Int]
assignPointsToClusters distanceMatrix medoids n = [findNearestMedoid i | i <- [0 .. n - 1]]
  where
    findNearestMedoid p = medoids !! minValIndex
      where
        distances = [distanceMatrix `atIndex` (p, i) | i <- medoids]
        minValIndex = minIndex (fromList distances)

-- | Calculate the average distance from a point to other points in the same cluster.
calculateAI :: Matrix Double -> Int -> Double
calculateAI subDistanceMatrix i
  | len <= 1 = 0
  | otherwise = sumElements (subDistanceMatrix ? [i]) / fromIntegral (len - 1)
  where
    len = rows subDistanceMatrix

-- | Calculate the minimum distance from a point to any point in other clusters.
calculateBI :: Matrix Double -> [Int] -> [Int] -> Int -> Double
calculateBI distanceMatrix clusters medoids i =
  let pointCluster = clusters !! i
      otherClusterMeans =
        [ meanDistance distanceMatrix otherClusterPoints i
          | medoid <- medoids,
            medoid /= pointCluster,
            let otherClusterPoints = getClusterPoints clusters medoid
        ]
   in if null otherClusterMeans
        then 1 / 0
        else minimum otherClusterMeans
        
-- | Calculate the mean distance from a point to all points in another cluster.
meanDistance :: Matrix Double -> [Int] -> Int -> Double
meanDistance distanceMatrix otherClusterPoints point
  | null otherClusterPoints = 1 / 0
  | otherwise = sumOtherDistances / fromIntegral (length otherClusterPoints)
  where
    sumOtherDistances = sumDistances distanceMatrix point otherClusterPoints

-- | Calculate silhouette score for each point in a cluster.
silhouetteForCluster :: Int -> Matrix Double -> [Int] -> [Int] -> Int -> [Double]
silhouetteForCluster n distanceMatrix clusters medoids medoidVal =
  let clusterPoints = getClusterPoints clusters medoidVal
      subDistanceMatrix = extractSubmatrix distanceMatrix clusterPoints
   in [silhouetteForPoint subDistanceMatrix clusterPoints j | j <- [0 .. length clusterPoints - 1]]
  where
    silhouetteForPoint subDistanceMatrix clusterPoints pointIndex =
      let a_i = calculateAI subDistanceMatrix pointIndex
          b_i = calculateBI distanceMatrix clusters medoids (clusterPoints !! pointIndex)
       in (b_i - a_i) / max a_i b_i

-- | Calculate the overall silhouette score for clustering.
silhouetteScore :: Int -> Matrix Double -> [Int] -> [Int] -> Double
silhouetteScore n distanceMatrix clusters medoids =
  let scores = silhouetteScoresRecursive medoids
   in mean scores
  where
    silhouetteScoresRecursive [] = []
    silhouetteScoresRecursive (medoid : rest) =
      silhouetteForCluster n distanceMatrix clusters medoids medoid ++ silhouetteScoresRecursive rest
    mean scoresArr = sum scoresArr / fromIntegral (length scoresArr)
    
-- | Execute the PAM algorithm to find the optimal clustering.
pamAlgorithm :: Int -> Matrix Double -> Int -> ([Int], [Int])
pamAlgorithm n distanceMatrix k = do
  iteratePAM n distanceMatrix clusters medoids 0
  where
    medoids = [0 .. k - 1] -- Choose medoids from 0 to k-1
    clusters = assignPointsToClusters distanceMatrix medoids n

-- | Iterate the PAM algorithm until clusters stabilize.
iteratePAM :: Int -> Matrix Double -> [Int] -> [Int] -> Int -> ([Int], [Int])
iteratePAM n distanceMatrix clusters medoids iteration =
  if clusters == newClusters
    then (newMedoids, newClusters)
    else iteratePAM n distanceMatrix newClusters newMedoids (iteration + 1)
  where
    newMedoids = findMedoid n distanceMatrix clusters medoids
    newClusters = assignPointsToClusters distanceMatrix newMedoids n

-- | Find the optimal number of clusters based on silhouette scores.
findOptimalClusters :: Int -> Matrix Double -> (Int, Double, [Int], [Int])
findOptimalClusters n distanceMatrix = findClusters 2 (-(1 / 0)) [] [] 0 2
  where
    maxIterationsWithoutImprovement = 5

    findClusters bestK bestScore bestMedoids bestClusters iterationsWithoutImprovement k
      | iterationsWithoutImprovement >= maxIterationsWithoutImprovement || k >= n =
          (bestK, bestScore, bestMedoids, bestClusters)
      | otherwise =
          let (medoids, clusters) = pamAlgorithm n distanceMatrix k
              silhouette = silhouetteScore n distanceMatrix clusters medoids
           in if silhouette > bestScore
                then findClusters k silhouette medoids clusters 0 (k + 1)
                else findClusters bestK bestScore bestMedoids bestClusters (iterationsWithoutImprovement + 1) (k + 1)

main :: IO ()
main = do
  -- Read matrix from the file
  ret <- readDistanceMatrix "dist_matrix(10000x10000).txt" ","
  let matrix = fromLists ret :: Matrix Double
  t_start <- getCurrentTime
  let n = rows matrix
  let (bestK, bestSilhouette, bestMedoids, bestClusters) = findOptimalClusters n matrix
  putStrLn $ "Optimal number of clusters: " ++ show bestK
  putStrLn $ "Best silhouette score: " ++ show bestSilhouette
  putStrLn $ "Best medoids: " ++ show bestMedoids
  t_end <- getCurrentTime
  putStrLn "PAM run time:"
  print $ diffUTCTime t_end t_start