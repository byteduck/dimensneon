import anndata
import numpy as np
from scipy.spatial import distance

# Based on the paper "Visualizing Data using t-SNE" by Laurens van der Maaten and Geoffrey Hinton
# https://www.jmlr.org/papers/volume9/vandermaaten08a/vandermaaten08a.pdf
# The following file contains an implementation of t-SNE. See the `tsne` function for main functionality.


def tsne(data: anndata.AnnData, perplexity=30.0, iterations=1000, animate=False):
    """
    Runs t-SNE on the given AnnData object.
    The t-SNE data will be embedded into the object in the frame "X_tsne".
    :param data: The AnnData object to run t-SNE on. PCA must have already been run on this object.
    :param perplexity: The perplexity to use for the simulation.
    :param iterations: The number of iterations to run.
    :param animate: Whether to store the datapoints as the simulation iterates for visualization purposes.
    """
    # The "X_pca" frame in our anndata object has our PCA data.
    X = data.obsm['X_pca']
    (n, dims) = X.shape

    # First, calculate the probability matrix for our data.
    # (That is, how likely it is for each point to interact with every other point)
    P = calculate_probability_matrix(X, perplexity)
    # As suggested by the paper, use "early exaggeration" for the first 1/4 iterations
    early_P = P * 4

    # Create an output array for the resultant 2-D points for our t-SNE plot. Start with random points
    Y = np.random.randn(n, 2)
    Ys = []
    momentum_Y = np.zeros((n, 2))
    learning_rates = np.ones((n, 2))
    LEARNING_RATE = 200

    for iter in range(iterations):
        if iter % 10 == 0:
            print(f"Running simulation (iteration {iter}/{iterations})...", end='\r')

        # If we're animating, we'll want to save the datapoints' positions so we can return
        # them. This was used for making a GIF for my presentation
        Ys.append(Y.copy())

        # First, calculate q_{i,j} matrix. (Equation 4)
        dists_inv = 1 / (1 + squared_distances(Y))
        dists_inv[range(n), range(n)] = 0.0 # Every point should have a distance of zero to itself
        Q = dists_inv / np.sum(dists_inv)

        # Now, it's gradient descent time!
        # Calculate the derivative of the Kullback-Leibler divergence between the probabilities (Equation 5)
        # (This is essentially just the change in momentum for each point)
        P_minus_Q = (early_P if iter < int(iterations * 0.25) else P) - Q
        delta_momentum = np.zeros((n, 2))
        for i in range(n):
            PQ_times_dist = np.tile(P_minus_Q[:, i] * dists_inv[:, i], (2, 1)).transpose()
            delta_momentum[i] = np.sum(PQ_times_dist * (Y[i] - Y), 0)

        # Where the momentum has increased, increase the learning rate. Where it's decreased, decrease it.
        # This isn't the exact method described in the paper, but the paper it referenced was a little out of my depth.
        # This seems to work as a nice approximation.
        learning_rates[(delta_momentum > 0.) != (momentum_Y > 0.)] += 0.1
        learning_rates[(delta_momentum > 0.) == (momentum_Y > 0.)] *= 0.9

        # Update the momentum of each point with our deltas in momentum according to the equation on pg. 4.
        alpha = 0.5 if iter < 250 else 0.8 # Alpha term of 0.5 for first 250 iterations, then 0.8 after
        momentum_Y = alpha * momentum_Y - LEARNING_RATE * (learning_rates * delta_momentum)
        Y += momentum_Y

    print(f"Running simulation (iteration {iterations}/{iterations})... Done!")

    # Finally, store the output in the "X_tsne" frame so that scanpy can plot it
    data.obsm['X_tsne'] = Y

    print("🎵 Now, these points of data make a beautiful line 🎵")
    print("🎵 And we're out of beta, we're releasing on time! 🎵")
    print("Finished!")

    if animate:
        return Ys


def squared_distances(X = np.array([[]])):
    """
    Calculates the squared Euclidean distance from every point in X to every other point in X.
    :param X: A numpy array of datapoints.
    :return: A numpy array of each datapoint's squared Euclidean distance to every other datapoint.
    """

    """
    # I realized that scipy has a faster implementation of this, which I switched to. The original code is below.
    # sqrt(Σ(a-b)^2) = sqrt(Σ(a^2 - 2ab + b^2)).

    # This is our a^2 & b^2 terms
    asquare = np.sum(np.square(X), 1)

    # This is our -2ab term
    twoab = -2 * np.dot(X, X.transpose())

    # This is Σ(a^2 - 2ab)
    diff = np.add(asquare, twoab).transpose()

    # Then finally, Σ((a^2-2ab) + b^2)
    squared_dists = np.add(diff, asquare)
    """

    return distance.squareform(distance.pdist(X, "sqeuclidean"))


def perplexity_probs(squared_dists, variance=1.0):
    """
    Calculates the Perplexity (normalized p_{i,j}) and normalized probability for a row of squared distances.
    :param squared_dists: The squared distances of every datapoint to the other.
    :param variance: The variance to use.
    :return: The perplexity of the row and normalized probabilities.
    """

    # Calculate the probability that each point chooses any other point as its neighbor, or p_{i,j}.
    # In other words, the probability based on a gaussian distribution at p_i at p_j with the given variance
    probs = np.exp(-squared_dists.copy() * variance)
    sum_probs = np.sum(probs) + 1e-15 # Avoid nonzero values

    # Normalize probs for perplexity
    norm_probs = probs / sum_probs

    # Calculate Perplexity of the row. We want to match this to the perplexity specified by the user
    perplexity = np.log(sum_probs) + variance * np.sum(squared_dists * probs) / sum_probs

    return perplexity, norm_probs


def calculate_probability_matrix(X = np.array([[]]), perplexity = 30.0):
    """
    Calculates the probability matrix for the given data.
    :param X: The datapoints to work on.
    :param perplexity: The target perplexity.
    :return: The probability matrix for the given data and perplexity.
    """
    (n, dims) = X.shape
    log_target_perplexity = np.log(perplexity)
    prob_matrix = np.zeros((n, n))

    # First, calculate distances from each point to the others
    print("Calculating squared Euclidean distances... ", end="")
    squared_dists = squared_distances(X)
    print("Done!")

    # Now, we want to hone in our variance for each row so that our perplexities match
    row_variance = np.ones((n, 1))
    for i in range(n):
        if i % 100 == 0:
            print(f"Calculating joint probabilities ({i}/{n})...", end='\r')

        # We want to take out X[i][i] from the row, since it's always 0
        col_indices = np.concatenate((np.r_[0:i], np.r_[i+1:n]))
        row = squared_dists[i, col_indices]

        # Hone in our variance so that we get a perplexity close to our target perplexity (arbitrary tolerance)
        # This is done with a binary search as suggested by the paper.
        TOLERANCE = 1e-4
        MAX_TRIES = 30
        max_var, min_var = (np.inf, -np.inf)
        tries = 0
        row_probs = np.zeros((n - 1))
        while tries < MAX_TRIES:
            # Calculate perplexity and probabilities for the row
            row_perplexity, row_probs = perplexity_probs(row, row_variance[i])
            perplexity_diff = row_perplexity - log_target_perplexity
            tries += 1

            # If the perplexity difference is within our tolerance, we're done!
            if np.abs(perplexity_diff) <= TOLERANCE:
                break

            # Otherwise, increase / decrease the variance as needed
            if perplexity_diff > 0:
                min_var = row_variance[i].copy()
                if max_var == np.inf or max_var == -np.inf:
                    row_variance[i] *= 2.0
                else:
                    row_variance[i] = (row_variance[i] + max_var) / 2.0
            else:
                max_var = row_variance[i].copy()
                if min_var == np.inf or min_var == -np.inf:
                    row_variance[i] /= 2.0
                else:
                    row_variance[i] = (row_variance[i] + min_var) / 2.0
        prob_matrix[i, col_indices] = row_probs
    print(f"Calculating joint probabilities ({n}/{n})... Done!")

    # Finally, normalize the probability matrix
    prob_matrix = prob_matrix + prob_matrix.transpose()
    prob_matrix = prob_matrix / np.sum(prob_matrix)

    return prob_matrix