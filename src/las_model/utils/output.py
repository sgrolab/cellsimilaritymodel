from las_model.utils.config import PROJECT_DIR
import pickle
import json 

def save_experiment(experiment_name, data, metadata, base_dir):
    """
    Saves experiment data + metadata into base_dir/experiment_name/
    """
    exp_dir = base_dir / experiment_name
    exp_dir.mkdir(parents=True, exist_ok=True)

    # Save the pickle
    pickle_path = exp_dir / f'{experiment_name}.pickle'
    with open(pickle_path, 'wb') as f:
        pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)

    # Save the metadata
    meta_path = exp_dir / 'metadata.json'
    with open(meta_path, 'w') as f:
        json.dump(metadata, f, indent=2)

    return exp_dir