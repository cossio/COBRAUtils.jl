import COBRA
import Utils # Pkg.clone("https://github.com/cossio/Utils.jl")

export readCobraModel, prepareXmlModel


"""
Read model in MATLAB format.
"""
function readCobraModel(path::String, matrix::String = "S", model::String = "model")
    COBRA.loadModel(path, matrix, model)
end


"""
Converts SBML (.xml) model to MATLAB (.mat) format.
"""
function prepareXmlModel(xmlpath::String)
    Utils.gunzip(xmlpath)
    matpath = Utils.removeExtension(xmlpath) * ".mat"
    @assert isfile(xmlpath) || isfile(matpath)
    if !isfile(matpath)
        pkgdir = Pkg.Dir.path() * "/MaxEntChemostat2018/"
        run(`python $(pkgdir * "/scripts/save_as_mat.py") $xmlpath`)
    end
    return matpath
end

