import glob
import pathlib
from pprint import pprint
import matplotlib.pyplot as plt
import numpy as np

import bluepyefe
import bluepyefe.extract
from bluepyefe.cell import Cell

files_metadata = {"219.2-E2000": {"IDRest": []}, "219.2-F1000": {"IDRest": []}, "219.2-F3001": {"IDRest": []},
            "219.2-G1000": {"IDRest": []}, "464.2-B2000": {"IDRest": []},
            "464.2-C1000": {"IDRest": []}}

# files_metadata = {"219.2-E2000": {"IDRest": []}}

for cell in files_metadata.keys():
    for file in glob.glob(f"./exp_data/{cell}.abf"):
        files_metadata[cell]["IDRest"].append({
            "filepath": file,
            "i_unit": "pA",
            "t_unit": "s",
            "v_unit": "mV",
            "ljp": 0.
        })

print(files_metadata)

cells = bluepyefe.extract.read_recordings(files_metadata=files_metadata)

pprint(vars(cells[0].recordings[0]))

interesting_efeatures = [
    'peak_time',
    'time_to_first_spike',
    'ISI_values',
    'irregularity_index',
    'adaptation_index',
    'spike_count_stimint',
    'AP_height',
    'min_AHP_values',
    'depolarized_base',
    'AP_duration_half_width',
    'AP_duration',
    'AP_rise_time',
    'AP_fall_time',
    'voltage_base',
    'decay_time_constant_after_stim',
    'minimum_voltage'
]

interesting_amplitudes = [0, 100, 200, 300]

targets = []
for efeature in interesting_efeatures:
    for amplitude in interesting_amplitudes:
        target = {
            "efeature": efeature,
            "protocol": "IDRest",
            "amplitude": amplitude,
            "tolerance": 20.,
        }

        targets.append(target)

efel_settings = {
    'strict_stiminterval': True,
    'Threshold': -20
}


cells = bluepyefe.extract.extract_efeatures_at_targets(
        cells=cells,
        targets=targets,
        efel_settings=efel_settings
)

# cells = bluepyefe.extract.extract_efeatures(
#         output_directory='./Figures/extractionFeatures',
#         files_metadata=files_metadata,
#         efel_settings=efel_settings
# )

for cell in cells:
    print("\nCell " + cell.name, ": ")
    [pprint(i.efeatures) for i in cell.recordings]

step=1
for cellId in range(6):
    for i in range(0, len(cells[cellId].recordings), step):
        plt.plot(cells[cellId].recordings[i].t, cells[cellId].recordings[i].voltage, label=str(i))
    plt.legend()
    # plt.show()

for cellId in range(6):
    for i in range(0, len(cells[cellId].recordings), step):
        plt.plot(cells[cellId].recordings[i].t, cells[cellId].recordings[i].current, label=str(i))
    plt.legend()
    # plt.show()

for cell in cells:
    cell.plot_all_recordings(show=False, output_dir='figures/PlotSteps')

Rheobase = {"219.2-E2000": 0.06, "219.2-F1000": 0.07, "219.2-F3001": 0.04,
            "219.2-G1000": 0.06, "464.2-B2000": 0.06,
            "464.2-C1000": 0.07}

for cell in cells:
    cell.rheobase = Rheobase[cell.name]
    cell.compute_relative_amp()

protocols = bluepyefe.extract.group_efeatures(cells, targets)

efeatures, protocol_definitions, currents = bluepyefe.extract.create_feature_protocol_files(
    cells=cells,
    protocols=protocols,
    output_directory='./figures/extractionFeatures'
)

pprint(efeatures)
pprint(protocol_definitions)
pprint(currents)
