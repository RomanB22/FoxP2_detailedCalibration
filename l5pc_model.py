"""Run simple cell optimisation"""

"""
Copyright (c) 2016-2020, EPFL/Blue Brain Project

 This file is part of BluePyOpt <https://github.com/BlueBrain/BluePyOpt>

 This library is free software; you can redistribute it and/or modify it under
 the terms of the GNU Lesser General Public License version 3.0 as published
 by the Free Software Foundation.

 This library is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
 details.

 You should have received a copy of the GNU Lesser General Public License
 along with this library; if not, write to the Free Software Foundation, Inc.,
 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
"""
# pylint: disable=R0914

import os
import json

import bluepyopt.ephys as ephys

script_dir = os.path.dirname(__file__)
config_dir = os.path.join(script_dir, 'config_3Comp')

# TODO store definition dicts in json
# TODO rename 'score' into 'objective'
# TODO add functionality to read settings of every object from config format


def define_mechanisms(mechanismSelected=[]):
    """Define mechanisms"""
    mech_definitions={}
    mech_definitionsOrig = load_mechanisms()

    for sec, mechs in mech_definitionsOrig.items():
        if sec=='all':
            mech_definitions[sec]=mech_definitionsOrig[sec]
        else:
            mech_definitions[sec]=[i for i in mechanismSelected]

    return create_mechanisms(mech_definitions)


def load_mechanisms():

    return json.load(
        open(
            os.path.join(
                config_dir,
                'mechanismsFull.json')))


def create_mechanisms(mech_definitions):

    mechanisms = []
    for sectionlist, channels in mech_definitions.items():
        seclist_loc = ephys.locations.NrnSeclistLocation(
            sectionlist,
            seclist_name=sectionlist)
        for channel in channels:
            mechanisms.append(ephys.mechanisms.NrnMODMechanism(
                name='%s.%s' % (channel, sectionlist),
                mod_path=None,
                suffix=channel,
                locations=[seclist_loc],
                preloaded=True))

    return mechanisms


def define_parameters(mechanismSelected=[]):
    """Define parameters"""
    param_configsOrig = load_parameters()
    param_configs = []

    for paramConf in param_configsOrig:
        if 'mech' in paramConf:
            if paramConf['mech'] in mechanismSelected:
                param_configs.append(paramConf)
        else:
            param_configs.append(paramConf)

    return create_parameters(param_configs)


def load_parameters():
    return json.load(open(os.path.join(config_dir, 'parametersFull.json')))


def create_parameters(param_configs):
    parameters = []

    for param_config in param_configs:
        if 'value' in param_config:
            frozen = True
            value = param_config['value']
            bounds = None
        elif 'bounds' in param_config:
            frozen = False
            bounds = param_config['bounds']
            value = None
        else:
            raise Exception(
                'Parameter config has to have bounds or value: %s'
                % param_config)

        if param_config['type'] == 'global':
            parameters.append(
                ephys.parameters.NrnGlobalParameter(
                    name=param_config['param_name'],
                    param_name=param_config['param_name'],
                    frozen=frozen,
                    bounds=bounds,
                    value=value))
        elif param_config['type'] in ['section', 'range']:
            if param_config['dist_type'] == 'uniform':
                scaler = ephys.parameterscalers.NrnSegmentLinearScaler()
            elif param_config['dist_type'] == 'exp':
                scaler = ephys.parameterscalers.NrnSegmentSomaDistanceScaler(
                    distribution=param_config['dist'])
            seclist_loc = ephys.locations.NrnSeclistLocation(
                param_config['sectionlist'],
                seclist_name=param_config['sectionlist'])

            name = '%s.%s' % (param_config['param_name'],
                              param_config['sectionlist'])

            if param_config['type'] == 'section':
                parameters.append(
                    ephys.parameters.NrnSectionParameter(
                        name=name,
                        param_name=param_config['param_name'],
                        value_scaler=scaler,
                        value=value,
                        frozen=frozen,
                        bounds=bounds,
                        locations=[seclist_loc]))
            elif param_config['type'] == 'range':
                parameters.append(
                    ephys.parameters.NrnRangeParameter(
                        name=name,
                        param_name=param_config['param_name'],
                        value_scaler=scaler,
                        value=value,
                        frozen=frozen,
                        bounds=bounds,
                        locations=[seclist_loc]))
        else:
            raise Exception(
                'Param config type has to be global, section or range: %s' %
                param_config)

    return parameters


def define_morphology(do_replace_axon, morphoFile = 'morphology/threeCompartmental.swc'):
    """Define morphology"""

    return ephys.morphologies.NrnFileMorphology(
        os.path.join(
            script_dir,
            morphoFile),
        do_replace_axon=do_replace_axon)


def create(do_replace_axon=True):
    """Create cell model"""

    cell = ephys.models.CellModel(
        'l5pc',
        morph=define_morphology(do_replace_axon),
        mechs=define_mechanisms(),
        params=define_parameters())

    return cell