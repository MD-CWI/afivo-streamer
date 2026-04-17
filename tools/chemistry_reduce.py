#!/usr/bin/env python3

# Tool for reduction of chemical reaction sets, using the "directed relation
# graph" method.
#
# Reference: Lu and Law, "A directed relation graph method for mechanism
# reduction" 2005

import numpy as np
from collections import deque

# Absolute tolerance for treating stoichiometric coefficients as zero
ATOL = 1e-30


def load_data(simname):
    """Load all simulation data files."""
    rates = np.loadtxt(f"{simname}_rates.txt")
    times = rates[:, 0]
    rates = rates[:, 1:]  # shape: (n_times, n_reactions)

    with open(simname + '_species.txt', 'r') as f:
        species = [x.strip() for x in f.readlines() if x.strip()]
    with open(simname + '_reactions.txt', 'r') as f:
        reactions = [x.strip() for x in f.readlines() if x.strip()]

    stoich_matrix = np.loadtxt(f"{simname}_stoich_matrix.txt")
    # shape: (n_species, n_reactions)

    n_species = len(species)
    n_reactions = len(reactions)

    # Load amounts of species in the total simulation volume
    amounts = np.loadtxt(simname + '_amounts.txt')
    amounts = amounts[:, 1:] # Drop time column

    assert stoich_matrix.shape == (n_species, n_reactions), \
        f"Stoichiometric matrix shape {stoich_matrix.shape} != " \
        f"({n_species}, {n_reactions})"
    assert rates.shape[1] == n_reactions, \
        f"Rates columns {rates.shape[1]} != n_reactions {n_reactions}"

    print(f"Loaded mechanism: {n_species} species, {n_reactions} reactions, "
          f"{rates.shape[0]} time steps")

    return {
        "times": times,
        "rates": rates,
        "species": species,
        "reactions": reactions,
        "stoich_matrix": stoich_matrix,
        "amounts": amounts
    }


def compute_drg_coefficients(stoich_matrix, rates):
    """
    Compute the DRG pair-wise coupling coefficients r_AB.

    For species A and B, the DRG coefficient is defined as:

        r_AB = sum_i |nu_{A,i} * R_i| * delta_Bi
               ----------------------------------------
               sum_i |nu_{A,i} * R_i|

    where:
      - nu_{A,i} is the stoichiometric coefficient of species A in reaction i
      - R_i is the rate of reaction i
      - delta_Bi = 1 if species B participates in reaction i, 0 otherwise

    We compute this for each time step and take the maximum over time.

    Parameters
    ----------
    stoich_matrix : array, shape (n_species, n_reactions)
        Net stoichiometric matrix.
    rates : array, shape (n_times, n_reactions)
        Reaction rates at each time step.

    Returns
    -------
    r_AB : array, shape (n_species, n_species)
        DRG coupling coefficients (maximum over all time steps).
    """
    participation = (np.abs(stoich_matrix) > ATOL).astype(float)

    # contrib[t, A, i] = |nu_{A,i} * R_i(t)|
    # shape: (n_times, n_species, n_reactions)
    contrib = np.abs(stoich_matrix[np.newaxis, :, :] * rates[:, np.newaxis, :])

    # Denominator: sum over reactions for each (time, species)
    denom = np.maximum(contrib.sum(axis=2), ATOL)  # (n_times, n_species)

    # numerator[t, A, B] = sum_i contrib[t, A, i] * participation[B, i]
    numerator = np.einsum('tsi,bi->tsb', contrib, participation)

    # r_AB = max over time of numerator / denom
    r_AB = (numerator / denom[:, :, np.newaxis]).max(axis=0)

    return r_AB


def bfs_reachable(r_AB, target_indices, threshold):
    """
    Find all species reachable from target species via DRG graph search.

    A species B is reachable if there exists a path from any target species
    to B such that every edge has r_AB >= threshold.

    Parameters
    ----------
    r_AB : array, shape (n_species, n_species)
        DRG coupling coefficients.
    target_indices : list of int
        Indices of target species.
    threshold : float
        Minimum coupling coefficient for an edge.

    Returns
    -------
    kept : set of int
        Indices of all reachable (kept) species.
    """
    n_species = r_AB.shape[0]
    kept = set(target_indices)
    queue = deque(target_indices)

    while queue:
        A = queue.popleft()
        for B in range(n_species):
            if B not in kept and r_AB[A, B] >= threshold:
                kept.add(B)
                queue.append(B)

    return kept


def kept_reactions_from_species(stoich_matrix, kept_species_set):
    """
    Determine which reactions to keep based on the kept species.

    A reaction is kept if ALL species involved in it are in the kept set.

    Parameters
    ----------
    stoich_matrix : array, shape (n_species, n_reactions)
        Net stoichiometric matrix.
    kept_species_set : set of int
        Indices of kept species.

    Returns
    -------
    kept_reaction_idx : list of int
        Indices of kept reactions.
    """
    kept = []
    for j in range(stoich_matrix.shape[1]):
        involved = np.flatnonzero(np.abs(stoich_matrix[:, j]) > ATOL)
        if all(s in kept_species_set for s in involved):
            kept.append(j)
    return kept


def drg_reduce(data, target_species, threshold=0.01, verbose=False):
    """
    Apply the DRG method to reduce the mechanism.

    Parameters
    ----------
    data : dict
        Data loaded by load_data().
    target_species : list of str
        Species that must be retained
    threshold : float
        DRG threshold epsilon. Species with coupling coefficient below this
        value to all target/retained species are removed.

    Returns
    -------
    result : dict
        Dictionary with keys:
        - kept_species : list of str
        - kept_reactions : list of str
        - kept_species_idx : list of int
        - kept_reactions_idx : list of int
        - r_AB : array, the DRG coefficient matrix
    """
    species = data["species"]
    reactions = data["reactions"]
    stoich_matrix = data["stoich_matrix"]
    rates = data["rates"]

    n_species = len(species)
    n_reactions = len(reactions)

    # Map species names to indices
    species_idx = {name: i for i, name in enumerate(species)}

    # Validate target species
    target_indices = []
    for sp in target_species:
        if sp not in species_idx:
            print(f"Warning: target species '{sp}' not found")
        else:
            target_indices.append(species_idx[sp])

    if not target_indices:
        raise ValueError("No valid target species provided!")

    r_AB = compute_drg_coefficients(stoich_matrix, rates)

    kept_species_set = bfs_reachable(r_AB, target_indices, threshold)
    kept_species_idx = sorted(kept_species_set)
    kept_species = [species[i] for i in kept_species_idx]

    kept_reactions_idx = kept_reactions_from_species(stoich_matrix,
                                                     kept_species_set)
    kept_reactions = [reactions[i] for i in kept_reactions_idx]

    if verbose:
        print(f"\nDRG Reduction Results (threshold={threshold}):")
        print(f"  Species: {n_species} -> {len(kept_species)} "
              f"(removed {n_species - len(kept_species)})")
        print(f"  Reactions: {n_reactions} -> {len(kept_reactions)} "
              f"(removed {n_reactions - len(kept_reactions)})")

    return {
        "kept_species": kept_species,
        "kept_reactions": kept_reactions,
        "kept_species_idx": kept_species_idx,
        "kept_reactions_idx": kept_reactions_idx,
        "r_AB": r_AB,
    }


def sweep_thresholds(data, target_species, thresholds=None):
    """
    Sweep over multiple DRG thresholds and report reduction for each.

    Parameters
    ----------
    data : dict
        Data loaded by load_data().
    target_species : list of str
        Species that must be retained.
    thresholds : list of float, optional
        Thresholds to test. Defaults to a logarithmic range.

    Returns
    -------
    results : list of dict
        Results for each threshold.
    """
    if thresholds is None:
        thresholds = np.logspace(-3, 0, 30)

    species = data["species"]
    stoich_matrix = data["stoich_matrix"]
    n_species = len(species)
    n_reactions = len(data["reactions"])

    species_idx = {name: i for i, name in enumerate(species)}
    target_indices = [species_idx[sp] for sp in target_species
                      if sp in species_idx]

    print("Computing DRG coefficients...")
    r_AB = compute_drg_coefficients(stoich_matrix, data["rates"])

    results = []
    print(f"\n{'Threshold':>10} | {'Species':>10} | {'Reactions':>10} | "
          f"{'% Species':>10} | {'% Reactions':>12}")
    print("-" * 65)

    for eps in thresholds:
        kept_species_set = bfs_reachable(r_AB, target_indices, eps)
        kept_species_idx = sorted(kept_species_set)
        kept_reactions_idx = kept_reactions_from_species(stoich_matrix,
                                                         kept_species_set)

        pct_sp = 100.0 * len(kept_species_idx) / n_species
        pct_rx = 100.0 * len(kept_reactions_idx) / n_reactions

        print(f"{eps:>10.4f} | {len(kept_species_idx):>10} | "
              f"{len(kept_reactions_idx):>10} | {pct_sp:>9.1f}% | "
              f"{pct_rx:>10.1f}%")

        results.append({
            "threshold": eps,
            "n_species": len(kept_species_idx),
            "n_reactions": len(kept_reactions_idx),
            "kept_species_idx": kept_species_idx,
            "kept_reactions_idx": kept_reactions_idx,
        })

    return results


def detect_qss_species(data, threshold=0.1, min_fraction=0.9,
                       exclude_species=None, verbose=False):
    """
    Detect species satisfying the quasi-steady-state approximation.

    A species is flagged as QSS if, for at least `min_fraction` of all
    time steps, its net production/consumption imbalance is below
    `threshold`.

    Parameters
    ----------
    data : dict
        Data loaded by load_data().
    threshold : float
        Maximum allowed relative imbalance delta_A for QSS.
    min_fraction : float
        Fraction of time steps that must satisfy the criterion
        (avoids flagging species that are only transiently in steady state).
    exclude_species : list of str, optional
        Species that should never be marked as QSS (e.g., targets,
        major species, electrons).

    Returns
    -------
    qss_species : list of str
        Species names identified as quasi-steady-state.
    delta : array, shape (n_times, n_species)
        The imbalance metric at each time step for each species.
    """
    if exclude_species is None:
        exclude_species = []

    stoich = data["stoich_matrix"]   # (n_species, n_reactions)
    rates = data["rates"]            # (n_times, n_reactions)
    species = data["species"]

    n_times, n_reactions = rates.shape
    n_species = len(species)

    # Contribution of each reaction to each species: nu_{A,i} * R_i(t)
    # shape: (n_times, n_species, n_reactions)
    omega = stoich[np.newaxis, :, :] * rates[:, np.newaxis, :]

    # Production: positive contributions; Consumption: negative contributions
    production = np.where(omega > 0, omega, 0.0).sum(axis=2)   # (n_times, n_species)
    consumption = np.where(omega < 0, -omega, 0.0).sum(axis=2) # (n_times, n_species)

    # Net rate and imbalance
    net = np.abs(production - consumption)
    scale = np.maximum(production, consumption)
    scale = np.maximum(scale, ATOL)  # avoid division by zero

    delta = net / scale  # (n_times, n_species)

    # A species is QSS if delta < threshold for >= min_fraction of time steps
    fraction_steady = (delta < threshold).sum(axis=0) / n_times  # (n_species,)

    qss_species = []

    if verbose:
        print(f"\n{'Species':>20} | {'Frac steady':>12} | {'Mean delta':>12} | "
              f"{'QSS?':>5}")
        print("-" * 60)

    for i, sp in enumerate(species):
        is_qss = (fraction_steady[i] >= min_fraction
                  and sp not in exclude_species)
        if is_qss:
            qss_species.append(sp)

        # Print species that are close to QSS for diagnostics
        if verbose and (fraction_steady[i] > 0.5 or is_qss):
            print(f"{sp:>20} | {fraction_steady[i]:>11.3f} | "
                  f"{delta[:, i].mean():>12.4e} | "
                  f"{'YES' if is_qss else 'no':>5}")

    if verbose:
        print(f"\nIdentified {len(qss_species)} QSS species (threshold="
              f"{threshold}, min_fraction={min_fraction})")

    return qss_species, delta


def drg_with_qss(data, target_species, threshold=0.1, qss_threshold=0.1,
                 qss_min_fraction=0.9, verbose=False):
    """
    Two-stage reduction:
    1. DRG to remove unimportant species
    2. QSS detection among the kept species to identify algebraic species
    """
    # Stage 1: DRG graph-based removal
    drg_result = drg_reduce(data, target_species, threshold)
    kept_species = drg_result["kept_species"]
    kept_reactions = drg_result["kept_reactions"]
    r_AB = drg_result["r_AB"]

    # Stage 2: Among kept species, find QSS candidates
    qss_candidates, delta = detect_qss_species(
        data,
        threshold=qss_threshold,
        min_fraction=qss_min_fraction,
        exclude_species=target_species,  # never eliminate targets
        verbose=verbose
    )

    # Only consider QSS species that survived DRG
    qss_in_kept = [sp for sp in qss_candidates if sp in kept_species]

    # Kept species that are not qss
    ode_species = [sp for sp in kept_species if sp not in qss_in_kept]

    if verbose:
        print(f"  Original:    {len(data['species'])} species, "
              f"{len(data['reactions'])} reactions")
        print(f"  After DRG:   {len(kept_species)} species, "
              f"{len(kept_reactions)} reactions")
        print(f"  QSS species: {len(qss_in_kept)} species (removed by QSS)")
        print(f"  ODE species: {len(ode_species)} species (kept after QSS)")

    return {
        "kept_species": kept_species,
        "kept_reactions": kept_reactions,
        "qss_species": qss_in_kept,
        "ode_species": ode_species,
        "r_AB": r_AB,
        "delta": delta,
    }


def main():
    import argparse

    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="DRG mechanism reduction for plasma chemistry")

    p.add_argument('simname', type=str,
                   help='Base name of simulation files (e.g. "sim" for '
                        'sim_rates.txt, sim_species.txt, etc.)')
    p.add_argument("-t", "--target", nargs="+", required=True,
                   help="Target species that must be retained")
    p.add_argument("-e", "--threshold", type=float, default=0.1,
                   help="DRG threshold")
    p.add_argument("-s", "--sweep", action="store_true",
                   help="Sweep over multiple thresholds")
    p.add_argument("-q", "--qss", action="store_true",
                   help="Enable QSS detection after DRG reduction")
    p.add_argument("--qss-threshold", type=float, default=0.1,
                   help="QSS imbalance threshold")
    p.add_argument("--qss-min-fraction", type=float, default=0.9,
                   help="Min. fraction of time steps satisfying QSS criterion")
    p.add_argument("-v", "--verbose", action="store_true",
                   help="Enable verbose output")

    args = p.parse_args()

    data = load_data(args.simname)

    if args.sweep:
        sweep_thresholds(data, args.target)
    elif args.qss:
        result = drg_with_qss(data, args.target,
                              threshold=args.threshold,
                              qss_threshold=args.qss_threshold,
                              verbose=args.verbose)

        print(f"\nODE species ({len(result['ode_species'])}):")
        for sp in result["ode_species"]:
            print(f"  {sp}")

        print(f"\nQSS species ({len(result['qss_species'])}):")
        for sp in result["qss_species"]:
            print(f"  {sp}")

        removed_species = sorted(
            set(data["species"]) - set(result["kept_species"]))
        print(f"\nRemoved species ({len(removed_species)}):")
        for sp in removed_species:
            print(f"  {sp}")

        print(f"\nKept reactions ({len(result['kept_reactions'])}):")
        for rx in result["kept_reactions"]:
            print(f"  {rx}")

        removed_reactions = sorted(
            set(data["reactions"]) - set(result["kept_reactions"]))
        print(f"\nRemoved reactions ({len(removed_reactions)}):")
        for rx in removed_reactions:
            print(f"  {rx}")
    else:
        result = drg_reduce(data, args.target, args.threshold,
                            verbose=args.verbose)
        kept_species = result["kept_species"]
        kept_reactions = result["kept_reactions"]

        print(f"\nKept species ({len(kept_species)}):")
        for sp in kept_species:
            print(f"  {sp}")

        removed_species = sorted(set(data["species"]) - set(kept_species))
        print(f"\nRemoved species ({len(removed_species)}):")
        for sp in removed_species:
            print(f"  {sp}")

        print(f"\nKept reactions ({len(kept_reactions)}):")
        for rx in kept_reactions:
            print(f"  {rx}")

        removed_reactions = sorted(
            set(data["reactions"]) - set(kept_reactions))
        print(f"\nRemoved reactions ({len(removed_reactions)}):")
        for rx in removed_reactions:
            print(f"  {rx}")


if __name__ == "__main__":
    main()
