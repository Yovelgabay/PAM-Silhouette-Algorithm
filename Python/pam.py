import numpy as np
import time

# Reads the distance matrix from a text file.
def read_distance_matrix(filepath, delimiter):
    return np.loadtxt(filepath, delimiter=delimiter)


# Assigns each point to the nearest medoid.
def assign_points_to_clusters(distance_matrix, medoids):
    n = distance_matrix.shape[0]
    clusters = np.empty(n, dtype=int)
    for p in range(n):
        distances = distance_matrix[p, medoids]
        min_val_index = np.argmin(distances)
        clusters[p] = medoids[min_val_index]
    return clusters


# Finds the optimal medoid for each cluster by minimizing the sum of distances within the cluster.
def find_medoid(distance_matrix, clusters, medoids):
    for m in range(len(medoids)):
        min_sum = np.inf
        min_point = -1
        points = np.where(clusters == medoids[m])[0]  # Get all points assigned to this medoid
        cluster_distances = distance_matrix[points][:, points]
        for i in range(len(cluster_distances)):
            sum_distances = np.sum(cluster_distances[i])
            if sum_distances < min_sum:
                min_sum = sum_distances
                min_point = points[i]
        medoids[m] = min_point
    return medoids


# Runs the PAM (Partitioning Around Medoids) algorithm to find medoids and cluster assignments.
def pam_algorithm(distance_matrix, k):
    n = distance_matrix.shape[0]

    # Initialize medoids from 0 to k-1
    medoids = np.arange(k)
    clusters = assign_points_to_clusters(distance_matrix, np.copy(medoids))
    iteration = 0
    while True:
        iteration += 1
        new_medoids = find_medoid(distance_matrix, np.copy(clusters), np.copy(medoids))
        new_clusters = assign_points_to_clusters(distance_matrix, new_medoids)
        # Check for convergence
        if np.array_equal(medoids, new_medoids):
            break
        medoids = new_medoids
        clusters = new_clusters
    return medoids, clusters


# Computes the silhouette score for the current clustering.
def silhouette_score(distance_matrix, clusters, medoids):
    n = distance_matrix.shape[0]
    silhouette_scores = np.zeros(n)

    for cluster_index in range(len(medoids)):
        cluster_points = np.where(clusters == medoids[cluster_index])[0]
        cluster_distances = distance_matrix[cluster_points][:, cluster_points]
        for i in range(len(cluster_points)):
            point = cluster_points[i]
            if len(cluster_points) > 1:
                a_i = np.sum(cluster_distances[i]) / (len(cluster_points) - 1)
            else:
                a_i = 0
            b_i = np.inf
            for other_cluster_index in range(len(medoids)):
                if other_cluster_index != cluster_index:
                    other_cluster_points = np.where(clusters == medoids[other_cluster_index])[0]
                    if len(other_cluster_points) > 0:
                        other_distances = distance_matrix[point][other_cluster_points]
                        mean_dist = np.mean(other_distances)
                        b_i = min(b_i, mean_dist)
            silhouette_scores[point] = (b_i - a_i) / max(a_i, b_i)

    return np.nanmean(silhouette_scores)


# Finds the optimal number of clusters by running PAM algorithm with varying number of clusters.
def find_optimal_clusters(distance_matrix):
    n = distance_matrix.shape[0]
    best_k = 2
    best_score = float('-inf')
    best_medoids = None
    best_clusters = None
    max_iterations_without_improvement = 5
    iterations_without_improvement = 0

    k = 2
    while iterations_without_improvement < max_iterations_without_improvement:
        if k >= n:
            return best_k, best_score, best_medoids, best_clusters
        medoids, clusters = pam_algorithm(distance_matrix, k)
        silhouette = silhouette_score(distance_matrix, clusters, medoids)
        print(f"Number of clusters: {k}, Silhouette Score: {silhouette}")
        if silhouette > best_score:
            best_score = silhouette
            best_k = k
            best_medoids = medoids
            best_clusters = clusters
            iterations_without_improvement = 0
        else:
            iterations_without_improvement += 1
        k += 1

    return best_k, best_score, best_medoids, best_clusters


# Custom PAM implementation
distance_matrix = read_distance_matrix(filepath='dist_matrix(10000x10000).txt', delimiter=',')
start_time = time.time()
best_k, best_score, best_medoids, best_clusters = find_optimal_clusters(distance_matrix)
end_time = time.time()
custom_pam_time = end_time - start_time
print(f"PAM Time: {custom_pam_time:.2f} seconds")
print(f"Optimal number of clusters: {best_k}")
print(f"Best Silhouette Score: {best_score:.4f}")
print(f"Best Medoids: {best_medoids}")

