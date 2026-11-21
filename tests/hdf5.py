"""The objects of an HDF5 file as stored, to compare two .metacool maps."""

import h5py
import numpy as np


def _plain(value):
    """A value read by h5py, comparable with ==: text as it is, and an array or a
    number as its dtype, shape and bytes (text of variable length in an array
    as a list)."""
    if isinstance(value, (str, bytes)):
        return value
    array = np.asarray(value)
    if array.dtype.kind == "O":
        return repr(array.tolist())
    return (array.dtype.str, array.shape, array.tobytes())


def walk(path):
    """Every group, dataset and attribute of an HDF5 file, as stored, by path: the
    order of the members of each group, the dtype, shape, chunks, filters and
    values of each dataset, the dtype and value of each attribute, and the soft
    links."""
    found = {}

    def attributes(node):
        return [(key, node.attrs.get_id(key).dtype.str, _plain(node.attrs[key])) for key in node.attrs]

    def visit(group, prefix):
        found[prefix or "/"] = ("group", list(group.keys()), attributes(group))
        for key in group:
            link = group.get(key, getlink=True)
            name = f"{prefix}/{key}"
            if isinstance(link, h5py.SoftLink):
                found[name] = ("soft link", link.path)
            elif isinstance(group[key], h5py.Group):
                visit(group[key], name)
            else:
                data = group[key]
                found[name] = (
                    data.dtype.str,
                    data.shape,
                    data.chunks,
                    data.compression,
                    data.compression_opts,
                    data.shuffle,
                    _plain(data[()]),
                    attributes(data),
                )

    with h5py.File(path, "r") as handle:
        visit(handle, "")
    return found
