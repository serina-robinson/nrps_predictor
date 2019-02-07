def is_invalid_sequence(sequences: Iterable[str]) -> bool:
    if not sequences:
        raise ValueError("No sequence provided")

    aa = {"A", "R", "N", "D", "B", "C", "E", "Q", "Z", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"}
    for seq in sequences:
        res = any(elem not in aa for elem in seq) # check if a non-amino acid value is included

    return res