begin
    include("ClickDetection.jl")
    using .ClickDetection
    as_dataframe([2, 3, 4], DateTime(2013))
end