**PAM + Silhouette Algorithm**

**Overview
**
This repository contains an implementation of the PAM (Partitioning Around Medoids) + Silhouette Algorithm for clustering, developed in Python and Haskell. The primary goal is to determine the optimal number of clusters in real time while following strict constraints on external library usage.

What is PAM & K-Medoids?

The k-medoids problem is a clustering problem similar to k-means. The name was coined by Leonard Kaufman and Peter J. Rousseeuw with their PAM (Partitioning Around Medoids) algorithm. Both the k-means and k-medoids algorithms are partitional (breaking the dataset up into groups) and attempt to minimize the distance between points labeled to be in a cluster and a point designated as the center of that cluster.

Unlike k-means, k-medoids chooses actual data points as centers (medoids or exemplars), which allows for greater interpretability of the cluster centers. K-medoids can be used with arbitrary dissimilarity measures, whereas k-means generally requires Euclidean distance for efficient solutions. Because k-medoids minimizes a sum of pairwise dissimilarities instead of a sum of squared Euclidean distances, it is more robust to noise and outliers than k-means.
