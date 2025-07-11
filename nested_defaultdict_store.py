from collections import defaultdict
import numpy as np

# GOAL: design a flexible & organized data storage structure for incoming CCEP/ERNA data, grouped by metadata parameters
# metadata parameters: amplitude, frequency, stim_channel

# define nested defaultdict tree: 
# amplitude -> frequency -> stim_channel -> list of data entries (example hierarchy)
metadata_keys = ['amplitude', 'frequency', 'stim_channel']

# define innermost level (list of data entries)
# recursive definition allows for infinite defaultdicts, any missing key automatically creates another defaultdict of the same type
def tree():
    return defaultdict(tree)

root = tree() # top-level storage (amplitude), leaf is list

# insert grouped data entries (call right after data is ready for storage (preprocessed))
def add_chunk(data, meta):
    node = root
    for k in metadata_keys:
        node = node[meta[k]]
    entry = {'data': data, 'trial_id': meta.get('trial_id')}    # sample metadata that comes with actual data
    node.setdefault('_chunks', []).append(entry)

# retrieve all data chunks that match full metadata query
def get_group(meta_query):
    node = root
    for k in metadata_keys:
        node = node.get(meta_query[k], {})
    return node.get('_chunks', [])

# enable partial querying
def get_partial(meta_query, node=None, level=0):
    if node is None:
        node = root

    results = []

    # if we are at leaf (list) level
    if level == len(metadata_keys):
        return node.get('_chunks', [])

    key = metadata_keys[level]

    if key in meta_query:
        # if current key is in query, descend into that subkey
        subnode = node.get(meta_query[key], {})
        results.extend(get_partial(meta_query, subnode, level + 1))
    else:
        # if current key is not in query, descend into all subkeys at this level
        for subkey in node:
            if subkey == '_chunks':
                continue        # skip if current key is not specified in query (ie a list/data) -> branch
            subnode = node[subkey]
            results.extend(get_partial(meta_query, subnode, level + 1))

    return results



# == TESTING ==
if __name__ == "__main__":
    dummy_data = np.random.randn(1500, 16)       # initialize 1500 samples, 16 channels timeseries data w/ random values
    dummy_meta = {'amplitude': 1.0, 'frequency': 130, 'stim_channel': 5, 'timestamp': 123456, 'trial_id': 1}
    add_chunk(dummy_data, dummy_meta)

    query_meta = {'amplitude': 1.0, 'frequency': 130, 'stim_channel': 5}
    results = get_group(query_meta)
    #print(results)

    # partial querying
    query = {'amplitude': 1.0}      # search for all lists with amplitude of 1
    results = get_partial(query)
    print("Found", len(results), "match(es) to query")
    print(results)

    query_2 = {'stim_channel': 5}
    results_2 = get_partial(query_2)
    print("Found", len(results_2), "match(es) to query")

