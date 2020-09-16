
def format_blast_params(params):
    return " ".join([f"-{key} {value}" for key, value in params.items()])
