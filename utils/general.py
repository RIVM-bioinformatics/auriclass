import argparse


def add_tag(tag, lines):
    """
    Prepend a tag to text explaining from which analysis the message came.
    """
    full_tag = f"[{tag}]"
    if lines:
        return "\n".join(
            [f"{full_tag} {line}" for line in lines.split("\n") if line != ""]
        )
    else:
        return f"{full_tag}"


def check_number_within_range(minimum=0, maximum=1):
    """
    Creates a function to check whether a numeric value is within a range, inclusive.

    The generated function can be used by the `type` parameter in argparse.ArgumentParser.
    See https://stackoverflow.com/a/53446178.

    Args:
        value: the numeric value to check.
        minimum: minimum of allowed range, inclusive.
        maximum: maximum of allowed range, inclusive.

    Returns:
        A function which takes a single argument and checks this against the range.

    Raises:
        argparse.ArgumentTypeError: if the value is outside the range.
        ValueError: if the value cannot be converted to float.
    """

    def generated_func_check_range(value: str) -> str:
        value_f = float(value)
        if (value_f < minimum) or (value_f > maximum):
            raise argparse.ArgumentTypeError(
                f"Supplied value {value} is not within expected range {minimum} to {maximum}."
            )
        return str(value)

    return generated_func_check_range
