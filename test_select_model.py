from src.openbind_model_selection.select_model import Pipeline, select_model, Constants

from pathlib import Path
import re
from rich import print as rprint

if __name__ == "__main__":
    path = Path('/dls/labxchem/data/lb42888/lb42888-6/processed/analysis/')
    dtag = 'Zika_NS5A-x5769'
    pandda_dir = path / 'pandda2' / 'panddas'
    pandda_dataset_dir = pandda_dir / 'processed_datasets' / dtag
    sts = []

    # Get the PanDDA 2 results
    sts.append(
            {
                'path': pandda_dataset_dir / 'modelled_structures' / f'{dtag}-pandda-model.pdb',
                'pipeline': Pipeline.PanDDA2,
                'pipeline_info': {
                    'pandda_path': pandda_dir,
                    'dtag': dtag
                }
            }
        )

    # Get pipedream results for each postrefine
    pipedream_dir = path / 'pipedream'
    pipedream_dataset_dir = pipedream_dir / dtag
    for st_model_dir in pipedream_dataset_dir.glob('postrefine*'):
        if st_model_dir.is_dir():
            match = re.match(Constants.POSTREFINE_LIG_CODE_REGEX, st_model_dir.name)
            sts.append(
                {
                    'path': st_model_dir / 'refine.pdb',
                    'pipeline': Pipeline.Pipedream,
                    'pipeline_info': {
                        'pipedream_path': pipedream_dir,
                        'rhofit_dir': pipedream_dataset_dir / f'rhofit-{re.match(Constants.POSTREFINE_LIG_CODE_REGEX, st_model_dir.name)[1]}'
                    }
                }
            )

    rprint(sts)

    result = select_model(
        sts
    )

    rprint('The result is: ')
    rprint(result)
    ...