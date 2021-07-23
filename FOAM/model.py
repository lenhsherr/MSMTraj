"""
Estimates markov models
"""
from typing import List, Generator, Union
from pathlib import Path

import deeptime.markov as markov
from mdtraj import compute_phi, compute_psi

import multiprocessing
import numpy as np
from glob import glob
import mdtraj as md
from simtk import openmm
from simtk.openmm import app
import os
import pyemma
import bz2

# n_procs = 16
# n_seeds = 50

# if you change these, you have to rm features.npy and recompute
#select_contacts = np.load('select_contacts.npy')

def calc_features(traj_path: Path, out_dir: Path) -> np.ndarray:
    """Featurize data, both with features used for the explore MSM and for the exploit goal.
    This function is used for a multiprocessing Pool.
    Parameters
    ----------
    traj_path : Path
        Path to an trajectory file.
    out_dir: Path
        Path to feature trajectory

    Returns
    -------
    features : 
        A 2D (framesxfeatures) numpy array of features
    """
    
    traj_name = traj_path.stem
    feat_path = out_dir.joinpath(f"{traj_name}.npy")

    if feat_path.exists()
        return np.load(str(feat_path), allow_pickle=True)
    else:
        traj = md.load(str(traj_path))
        features = md.compute_contacts(traj)[0]
        np.save(file=feat_path.open('wb'), arr=features)
    return features

def calc_msm_counts(features):
    """Do tICA and k-means clustering.
    Parameters
    ----------
    features : list
        List of numpy feature arrays
    Returns
    -------
    counts : numpy array
        Array of counts each microstate was seen
    dtrajs : list
        List of discretized trajectories
    """
    tica = pyemma.coordinates.tica(features, lag=10, dim=5, kinetic_map=True, commute_map=False)
    tica_projection = tica.get_output()
    if len(np.concatenate(tica_projection)) > 100:
        kmeans = pyemma.coordinates.cluster_kmeans(tica_projection, k=100, max_iter=1000)
    else: # for beginning of sampling
        kmeans = pyemma.coordinates.cluster_kmeans(tica_projection, k=10, max_iter=1000)
    dtrajs = kmeans.dtrajs
    counts = np.bincount(np.concatenate(dtrajs))
    return counts, dtrajs

def draw_seeds(goals, counts, dtrajs, n_seeds, explore_weight=0.5, truncate_prob=True, max_states=None):
    """Calculate microstate scores, draw new seed states and random frames from them.
    Parameters
    ----------
    goals : list
        List of trajectories scored with one chosen goal feature
    counts : numpy array
        Array of counts each microstate was seen
    dtrajs : list
        List of discretized trajectories
    n_seeds : int
        Number of new seeds to produce
    explore_weight : float
        Weight of the explore part in the FAST scoring function
    truncate_prob : bool
        Whether to set to 0 the draw probabilities of states below 
        a number (determined by heuristic copied from HTMD) of top scoring states
    max_states : None or int
        If truncate_prob, maximum number of states to retain non-zero draw probabilities
    Returns
    -------
    seeds : numpy array
        array of of [trajectory, frame] indexes of new seeds
    """
    indexes = []
    for i,traj in enumerate(goals):
        for j,frame in enumerate(traj):
            indexes.append([i,j])
    indexes = np.array(indexes)

    goals = np.concatenate(goals)
    dtrajs = np.concatenate(dtrajs)
    
    goal_scores = []
    for state in range(len(set(dtrajs))):
        goal = np.mean(goals[np.argwhere(dtrajs == state).flatten()])
        goal_scores.append(goal)
    # normalize to [0,1]
    goal_scores = (goal_scores - np.min(goal_scores)) / (np.max(goal_scores) - np.min(goal_scores))

    count_scores = 1/counts
    # normalize to [0,1] (if all counts the same, set to uniform)
    if np.max(count_scores) - np.min(count_scores) == 0:
        count_scores = np.ones(len(count_scores)) / len(count_scores)
    else:
        count_scores = (count_scores - np.min(count_scores)) / (np.max(count_scores) - np.min(count_scores))

    combo_scores = (1-explore_weight)*goal_scores + explore_weight*count_scores
    if truncate_prob: # using the cumsum method from htmd
        order = np.argsort(combo_scores)[::-1]
        combo_scores_sorted = combo_scores[order]
        truncate_bool = (n_seeds * combo_scores_sorted / np.cumsum(combo_scores_sorted)) < 1
        if not max_states:
            combo_scores[order[truncate_bool]] = 0
        else:
            if len(truncate_bool) - sum(truncate_bool) > max_states:
                combo_scores[order[max_states:]] = 0
            else:
                combo_scores[order[truncate_bool]] = 0
    # normalize to sum 1
    combo_scores = combo_scores / np.sum(combo_scores)
    
    # draw states
    states_drawn = np.random.choice(range(len(set(dtrajs))), n_seeds, replace=True, p=combo_scores)
    
    # randomly draw frames from those states
    seeds = []
    for state in states_drawn:
        indexes_ = indexes[np.argwhere(dtrajs == state).flatten()]
        seeds_ = np.random.choice(range(len(indexes_)), 1)
        seeds_ = indexes_[seeds_]
        seeds.append(seeds_)
    seeds = np.concatenate(seeds)

    return seeds
    
def make_seed(trajs: List[Path], index: Tuple[int, int], seed_num: int, seed_dir: Path):
    """Save OpenMM state.xml file of a new seed with random new velocities.
    This function is used for a multiprocessing Pool.
    Parameters
    ----------
    trajs : list
        List of trajectory paths
    index : numpy array
        [trajectory, frame] indexes of the seed
    i : int
        Number of the seed to write as i.xml
    """
    seed = md.load_frame(trajs[index[0]], frame=index[1])
    seed.save('seeds_pdb/%d.pdb' % i)
    box_vectors = seed.openmm_boxes(0)

    pdb = app.PDBFile('seeds_pdb/%d.pdb' % i)
    
    with open('system.xml', 'r') as infile:
        system = openmm.XmlSerializer.deserialize(infile.read())

    with open('integrator.xml', 'r') as infile:
        integrator = openmm.XmlSerializer.deserialize(infile.read())

    context = openmm.Context(system, integrator)

    context.setPeriodicBoxVectors(*box_vectors)
    context.setPositions(pdb.positions)
    context.setVelocitiesToTemperature(310)
    context.applyConstraints(integrator.getConstraintTolerance())
    context.applyVelocityConstraints(integrator.getConstraintTolerance())

    state = context.getState(getPositions=True, getVelocities=True, getEnergy=True, getForces=True)

    with bz2.open('SEEDS/%d.xml.bz2' % i, 'wt') as outfile:
        xml = openmm.XmlSerializer.serialize(state)
        outfile.write(xml)


def generator_length(g: Generator) -> int:
    return sum(1 for x in g)


def get_trajectories(traj_dir: Path, filter: str) -> Union[None, List[md.Trajectory]]: 
    traj_paths = traj_dir.glob(filter)
    if generator_length(traj_paths) == 0:
        return None
    return [md.load(str(x)) for x in traj_paths]






# def main():

#     traj_dir = '/home/rob/FOAM/scratch/trajectories'
#     feat_dir = '/home/rob/FOAM/scratch/features'

#     trajs = get_trajectories(traj_dir, filter='*.h5')
    
#     # features = pool.map(calc_features, trajs)
#     # # to speed up tica let's take random 1000 features
#     # choice = np.random.choice(range(features[0].shape[1]), 1000)
#     # counts, dtrajs = calc_msm_counts([x[:,choice] for x in features])
#     # # save dtrajs for checking progress
#     # np.save('current_dtrajs.npy', dtrajs)
#     # np.save('current_trajs_glob.npy', trajs)
#     # # explore only - no goal, just pass dtrajs to satisfy the argument
#     # goals = dtrajs

#     # # with open('current_goal.txt', 'w') as f:
#     # #     f.write('Current max. goal: %f' % np.max(goals))
#     # #     f.write('\n')

#     # # explore only - set explore weight to zero, 
#     # seeds = draw_seeds(goals, counts, dtrajs, n_seeds, explore_weight=0, truncate_prob=True)

#     # pool.starmap(make_seed, [(trajs, seeds[i], i) for i in range(len(seeds))])

#     # os.system('sudo cp -r SEEDS /home/server/server2/projects/17805/')

# if __name__ == "__main__":
#     # pool = multiprocessing.Pool(n_procs)
#     # while True:
#     #     main()
