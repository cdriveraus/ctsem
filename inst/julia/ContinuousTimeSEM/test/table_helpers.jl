# DataFrame conveniences for the test suite only.
#
# The package itself takes the parameter table as plain column vectors (see
# `src/r_interface.jl` for why DataFrames is no longer a dependency of it), but
# a `DataFrame` literal is still the clearest way to write a small table in a
# test. These helpers convert one, translating `missing` into the sentinels the
# column API expects.

using DataFrames

_tbl_column(df, name, default) =
    name in propertynames(df) ? getproperty(df, name) : fill(default, nrow(df))

_tbl_int(column) = [ismissing(v) ? 0 : Int(v) for v in column]
_tbl_float(column) = [ismissing(v) ? NaN : Float64(v) for v in column]
_tbl_string(column) = [ismissing(v) ? "" : String(v) for v in column]

"""
    ekf_from_data_frame(df, ti_effects=..., diffusion_state_indices=Int[];
                        continuous_time=true)

Test-only wrapper around `ContinuousTimeSEM.ekf_from_columns`, preserving the
call shape the tests were written against.
"""
function ekf_from_data_frame(df::DataFrame,
    ti_effects::DataFrame=DataFrame(parameter=Int[], predictor=Int[], coefficient=Int[]),
    diffusion_state_indices::AbstractVector=Int[]; continuous_time::Bool=true)
    return ContinuousTimeSEM.ekf_from_columns(
        df.matrix, df.row, df.col,
        _tbl_int(_tbl_column(df, :parnumber, missing)),
        _tbl_float(_tbl_column(df, :value, missing)),
        _tbl_string(_tbl_column(df, :transform, missing)),
        _tbl_string(_tbl_column(df, :predicttransform, missing)),
        _tbl_string(_tbl_column(df, :updatetransform, missing)),
        _tbl_string(_tbl_column(df, :tdtransform, missing)),
        ti_parameter=Int.(ti_effects.parameter),
        ti_predictor=Int.(ti_effects.predictor),
        ti_coefficient=Int.(ti_effects.coefficient),
        diffusion_state_indices=Int.(diffusion_state_indices),
        continuous_time=continuous_time)
end

"""Test-only wrapper: axis for a layout given as a `DataFrame`."""
retrieve_axes(df::DataFrame) = ContinuousTimeSEM.retrieve_axes(df.matrix, df.row, df.col)

"""Test-only wrapper: names and dims for a layout given as a `DataFrame`."""
retrieve_names_and_dims(df::DataFrame) =
    ContinuousTimeSEM.retrieve_names_and_dims(df.matrix, df.row, df.col)
