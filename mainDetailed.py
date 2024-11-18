import bluepyopt as bpopt
import bluepyopt.ephys as ephys
import matplotlib.pyplot as plt

plotMorpho=False
plotResponses=True
showResponses=True
verbose=False
offspringSize=50
maxGenerations=20

morphoFile = './morphology/threeCompartmental.swc' # threeCompartmental
# threeCompartmental characteristics
# Soma surface area: 1256.6370614359173 um**2
# Neuron surface area: 4932.3004661359755 um**2


workDir = './config_3Comp/'
# mechanismOriginal = ['NaTs', 'Nap', 'Kv3_1', 'K_T', 'K_P', 'Ih', 'Im', 'Ca_HVA', 'Ca_LVA', 'SK']
# mechanismSelected = ['Nafx', 'Nap', 'kdrin', 'K_T', 'K_P', 'Ih', 'Im', 'Ca_HVA', 'Ca_LVA', 'SK']
mechanismSelected = ['NaV', 'Kv3_1']

if plotMorpho:
    import neurom
    from neurom.view import plot_morph, plot_morph3d
    plot_morph3d(neurom.load_morphology(morphoFile)) #
    plt.show()

    from neurom import features
    import numpy as np
    m = neurom.load_morphology(morphoFile)
    SomaSurf = features.get('soma_surface_area', m)
    NeuritesSurf = np.sum(features.get('total_area_per_neurite', m))
    print(features.get('total_area_per_neurite', m),features.get('total_length_per_neurite', m))
    print("Soma surface area: " + str(SomaSurf) + " um^2", '\n', "Neuron surface area: " + str(SomaSurf+NeuritesSurf) + " um^2")
    quit()
morphology = ephys.morphologies.NrnFileMorphology(morphology_path=morphoFile, do_replace_axon=True)

import json
param_configs = json.load(open(workDir+'parameters.json'))
if verbose: print([param_config['param_name'] for param_config in param_configs])

import l5pc_model
parameters = l5pc_model.define_parameters(mechanismSelected=mechanismSelected)
if verbose: print('\n'.join('%s' % param for param in parameters))

mechanisms = l5pc_model.define_mechanisms(mechanismSelected=mechanismSelected)
if verbose: print('\n'.join('%s' % mech for mech in mechanisms))

l5pc_cell = ephys.models.CellModel('l5pc', morph=morphology, mechs=mechanisms, params=parameters)
if verbose: print(l5pc_cell)

param_names = [param.name for param in l5pc_cell.params.values() if not param.frozen]

proto_configs = json.load(open(workDir+'protocols.json'))
if verbose: print(proto_configs)

import l5pc_evaluator
fitness_protocols = l5pc_evaluator.define_protocols()
if verbose: print('\n'.join('%s' % protocol for protocol in fitness_protocols.values()))

feature_configs = json.load(open(workDir+'features.json'))
if verbose: print(feature_configs)

fitness_calculator = l5pc_evaluator.define_fitness_calculator(fitness_protocols)
if verbose: print(fitness_calculator)

sim = ephys.simulators.NrnSimulator()

evaluator = ephys.evaluators.CellEvaluator(
        cell_model=l5pc_cell,
        param_names=param_names,
        fitness_protocols=fitness_protocols,
        fitness_calculator=fitness_calculator,
        sim=sim)

def plot_responses(responses, filename='./figures/optResults/responses.png'):
    fig, axes = plt.subplots(len(responses), figsize=(10,8))
    for index, (resp_name, response) in enumerate(responses.items()):
        axes[index].plot(response['time'], response['voltage'], label=resp_name)
        axes[index].set_title(resp_name)
    fig.tight_layout()
    fig.savefig(filename)

opt = bpopt.deapext.optimisationsCMA.DEAPOptimisationCMA(
    evaluator=evaluator,
    offspring_size=offspringSize,
    selector_name="multi_objective") #single_objective multi_objective

final_pop, halloffame, log, hist = opt.run(max_ngen=maxGenerations, cp_filename='checkpoints/checkpoint.pkl')

for i in range(len(halloffame)):
    best_params = evaluator.param_dict(halloffame[i])
    print('Best %d solution: \n' % i + str(best_params))

    best_responses = evaluator.run_protocols(protocols=fitness_protocols.values(), param_values=best_params)
    if plotResponses:
        plot_responses(best_responses, filename='./figures/optResults/Best_responses_%i.png' % i)
        if showResponses: plt.show()
