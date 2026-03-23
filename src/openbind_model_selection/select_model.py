from typing import TypedDict
from pathlib import Path
from enum import Enum
import re
import pandas as pd
import gemmi

Pipeline = Enum('Pipeline', ['PanDDA2', 'Pipedream', ])

PanDDA2Info = TypedDict(
    'StructureModel', 
    {
        'pandda_path': str | Path,
        'dtag': str
    }
)

PipedreamInfo = TypedDict(
    'PipedreamInfo', 
    {
        'pipedream_path': str | Path,
        'rhofit_dir': str | Path
    }
)

StructureModel = TypedDict(
    'StructureModel', 
    {
        'path': str | Path,
        'pipeline': Pipeline,
        'pipeline_info': PanDDA2Info | PipedreamInfo 
        }
    )
class Constants:
    RHOFIT_HIT_LOG =  'Hit_corr.log'
    PANDDA_ANALYZE_EVENTS = 'panddas/analyses/pandda_analyse_events.csv'
    PANDDA_ANALYZE_EVENT_SCORE = 'z_mean'
    LIG_NAMES = ('LIG', 'XXX')
    POSTREFINE_LIG_CODE_REGEX = r'postrefine\-([\S]+)'


def has_ligand(st: StructureModel) -> bool:
    struct = gemmi.read_structure(str(st['path']))
    for model in struct:
        for chain in model:
            for residue in chain:
                if residue.name in Constants.LIG_NAMES:
                    return True

    return False
    ...


def get_rscc_rhofit(st: StructureModel) -> float:
    # Actually gets the best RSCC of any fit, not -strictly- the one rhofit chose
    with open(Path(st['pipeline_info']['rhofit_dir']) / Constants.RHOFIT_HIT_LOG, 'r') as f:
        data = f.read()
    matches = re.findall(r'^[\S]\s([\S]+)', data)

    return max([float(match) for match in matches])

def get_pandda_score(st: StructureModel) -> float:
    df = pd.read_csv(Path(st['pipeline_info']['pandda_path'] / Constants.PANDDA_ANALYZE_EVENTS))

    df_dtag = df[df['dtag'] == st['pipeline_info']['dtag']]
    highest_event_score = df_dtag[Constants.PANDDA_ANALYZE_EVENT_SCORE].max()
    return highest_event_score


def select_model(sts: list[StructureModel]) -> StructureModel | None:
    """
    Given a list of models, produced by various pipelines, select which to continue
    """
    # Consider only structures with built ligands
    sts_with_ligands = {
        pipeline: {x['path']: x for x in sts if (x['pipeline'] == pipeline) & (has_ligand(x))}
        for pipeline 
        in Pipeline
        
    }

    if Pipeline.PanDDA2 in sts_with_ligands & Pipeline.Pipedream in sts_with_ligands:

        # If there are pipedream options, use the one with the best RSCC
        if len(sts_with_ligands[Pipeline.Pipedream]) > 0:
            rsccs = {path: get_rscc_rhofit(st) for path, st in sts_with_ligands[Pipeline.Pipedream].items()}
            return sts_with_ligands[Pipeline.Pipedream][max(rsccs, key=lambda _path: rsccs[_path])]

        # Otherwise if there is a PanDDA option, use the one with the highest event score
        elif len(sts_with_ligands[Pipeline.PanDDA2]) > 0:
            scores = {path: get_pandda_score(st) for path, st in sts_with_ligands[Pipeline.PanDDA2].items()}
            return sts_with_ligands[Pipeline.PanDDA2][max(scores, key=lambda _path: scores[_path])]

        # If there are no options return None
        else:
            return None
        
    elif Pipeline.Pipedream in sts_with_ligands:
        # Select pipedream with best RSCC
        rsccs = {path: get_rscc_rhofit(st) for path, st in sts_with_ligands[Pipeline.Pipedream].items()}
        return sts_with_ligands[Pipeline.Pipedream][max(rsccs, key=lambda _path: rsccs[_path])]
        
    elif Pipeline.PanDDA2 in sts_with_ligands:
        # Select highest event score PanDDA
        scores = {path: get_pandda_score(st) for path, st in sts_with_ligands[Pipeline.PanDDA2].items()}
        return sts_with_ligands[Pipeline.PanDDA2][max(scores, key=lambda _path: scores[_path])]        

    else:
        # No valid ligand
        return None
        ...
    ...