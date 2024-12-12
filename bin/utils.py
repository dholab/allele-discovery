def removeSpecialCharacters(in_str, special_characters="*|: ", replace_character="_"):
    """remove specified special characters from input str"""

    return in_str.translate({ord(c): replace_character for c in special_characters})
