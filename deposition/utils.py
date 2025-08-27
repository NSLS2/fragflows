import string


def letter_generator():
    def get_letter(num):
        result = []
        while num > 0:
            num, remainder = divmod(num - 1, 26)
            result.append(string.ascii_uppercase[remainder])
            return "".join(reversed(result))

    for k in range(1, 15000):
        yield get_letter(k)

def flatten_dict(d, parent_key="", sep="."):
    items = []
    for k, v in d.items():
        new_key = f"{parent_key}{sep}{k}" if parent_key else k
        if isinstance(v, dict):
            items.extend(flatten_dict(v, new_key, sep=sep).items())
        else:
            items.append((new_key, v))
    return dict(items)