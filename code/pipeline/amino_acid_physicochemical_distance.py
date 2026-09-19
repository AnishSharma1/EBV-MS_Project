"""Published amino-acid physicochemical distances for exact 9-mer comparisons."""

from __future__ import annotations

import math
from typing import Mapping, Sequence


TCR_FACING_INDICES = (1, 2, 4, 6, 7)

# Table 2 of Grantham R. Science. 1974;185:862-864.
# DOI: 10.1126/science.185.4154.862
_GRANTHAM_ORDER = "SRLPTAVGIFYCHQNKDEMW"
_GRANTHAM_ROWS = """
0 110 145 74 58 99 124 56 142 155 144 112 89 68 46 121 65 80 135 177
110 0 102 103 71 112 96 125 97 97 77 180 29 43 86 26 96 54 91 101
145 102 0 98 92 96 32 138 5 22 36 198 99 113 153 107 172 138 15 61
74 103 98 0 38 27 68 42 95 114 110 169 77 76 91 103 108 93 87 147
58 71 92 38 0 58 69 59 89 103 92 149 47 42 65 78 85 65 81 128
99 112 96 27 58 0 64 60 94 113 112 195 86 91 111 106 126 107 84 148
124 96 32 68 69 64 0 109 29 50 55 192 84 96 133 97 152 121 21 88
56 125 138 42 59 60 109 0 135 153 147 159 98 87 80 127 94 98 127 184
142 97 5 95 89 94 29 135 0 21 33 198 94 109 149 102 168 134 10 61
155 97 22 114 103 113 50 153 21 0 22 205 100 116 158 102 177 140 28 40
144 77 36 110 92 112 55 147 33 22 0 194 83 99 143 85 160 122 36 37
112 180 198 169 149 195 192 159 198 205 194 0 174 154 139 202 154 170 196 215
89 29 99 77 47 86 84 98 94 100 83 174 0 24 68 32 81 40 87 115
68 43 113 76 42 91 96 87 109 116 99 154 24 0 46 53 61 29 101 130
46 86 153 91 65 111 133 80 149 158 143 139 68 46 0 94 23 42 142 174
121 26 107 103 78 106 97 127 102 102 85 202 32 53 94 0 101 56 95 110
65 96 172 108 85 126 152 94 168 177 160 154 81 61 23 101 0 45 160 181
80 54 138 93 65 107 121 98 134 140 122 170 40 29 42 56 45 0 126 152
135 91 15 87 81 84 21 127 10 28 36 196 87 101 142 95 160 126 0 67
177 101 61 147 128 148 88 184 61 40 37 215 115 130 174 110 181 152 67 0
"""
_GRANTHAM_MATRIX = {
    (left, right): value
    for left, row in zip(_GRANTHAM_ORDER, _GRANTHAM_ROWS.strip().splitlines())
    for right, value in zip(_GRANTHAM_ORDER, (int(item) for item in row.split()))
}

# Table 2 of Atchley WR et al. PNAS. 2005;102:6395-6400.
# DOI: 10.1073/pnas.0408677102
ATCHLEY_FACTORS: Mapping[str, tuple[float, float, float, float, float]] = {
    "A": (-0.591, -1.302, -0.733, 1.570, -0.146),
    "C": (-1.343, 0.465, -0.862, -1.020, -0.255),
    "D": (1.050, 0.302, -3.656, -0.259, -3.242),
    "E": (1.357, -1.453, 1.477, 0.113, -0.837),
    "F": (-1.006, -0.590, 1.891, -0.397, 0.412),
    "G": (-0.384, 1.652, 1.330, 1.045, 2.064),
    "H": (0.336, -0.417, -1.673, -1.474, -0.078),
    "I": (-1.239, -0.547, 2.131, 0.393, 0.816),
    "K": (1.831, -0.561, 0.533, -0.277, 1.648),
    "L": (-1.019, -0.987, -1.505, 1.266, -0.912),
    "M": (-0.663, -1.524, 2.219, -1.005, 1.212),
    "N": (0.945, 0.828, 1.299, -0.169, 0.933),
    "P": (0.189, 2.081, -1.628, 0.421, -1.392),
    "Q": (0.931, -0.179, -3.005, -0.503, -1.853),
    "R": (1.538, -0.055, 1.502, 0.440, 2.897),
    "S": (-0.228, 1.399, -4.760, 0.670, -2.647),
    "T": (-0.032, 0.326, 2.213, 0.908, 1.313),
    "V": (-1.337, -0.279, -0.544, 1.242, -1.262),
    "W": (-0.595, 0.009, 0.672, -2.128, -0.184),
    "Y": (0.260, 0.830, 3.097, -0.838, 1.512),
}


def _validate_residue(aa: str) -> None:
    if aa not in ATCHLEY_FACTORS:
        raise ValueError(f"unsupported amino acid: {aa}")


def _validate_cores(left_core: str, right_core: str) -> None:
    if len(left_core) != 9 or len(right_core) != 9:
        raise ValueError("comparison requires two exact nine-residue cores")
    for aa in left_core + right_core:
        _validate_residue(aa)


def grantham_residue_distance(left: str, right: str) -> int:
    """Return the integer distance published in Grantham Table 2."""
    _validate_residue(left)
    _validate_residue(right)
    return _GRANTHAM_MATRIX[(left, right)]


def atchley_residue_distance(left: str, right: str) -> float:
    """Return ordinary Euclidean distance across the five published factors."""
    _validate_residue(left)
    _validate_residue(right)
    return math.sqrt(sum(
        (left_value - right_value) ** 2
        for left_value, right_value in zip(ATCHLEY_FACTORS[left], ATCHLEY_FACTORS[right])
    ))


def _mean_position_distance(
    left_core: str,
    right_core: str,
    distance,
    positions: Sequence[int] = TCR_FACING_INDICES,
) -> float:
    _validate_cores(left_core, right_core)
    indices = tuple(int(index) for index in positions)
    if not indices or any(index < 0 or index >= 9 for index in indices):
        raise ValueError("comparison positions are out of range")
    return sum(distance(left_core[index], right_core[index]) for index in indices) / len(indices)


def tcr_face_grantham_mismatch(left_core: str, right_core: str) -> float:
    """Mean Grantham distance across P2, P3, P5, P7, and P8."""
    return _mean_position_distance(left_core, right_core, grantham_residue_distance)


def tcr_face_atchley_distance(left_core: str, right_core: str) -> float:
    """Mean five-factor Euclidean distance across P2, P3, P5, P7, and P8."""
    return _mean_position_distance(left_core, right_core, atchley_residue_distance)
