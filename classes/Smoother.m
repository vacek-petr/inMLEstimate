classdef Smoother
% SMOOTHER class for calling smoothing routines

properties
    name
    numberOfPreSmoothingIterations  % vector - possible to choose different number on each level
    numberOfPostSmoothingIterations % vector - possible to choose different number on each level
    dampingParameter % number - supported in jacobi, sor
end

methods
    function obj = Smoother(name,numberOfPreSmoothingIterations,numberOfPostSmoothingIterations,options)
    % SMOOTHER Construct an instance of class Smoother
    
    arguments
        name {mustBeMember(name,{'jacobi','gs'})}
        numberOfPreSmoothingIterations
        numberOfPostSmoothingIterations
        options.dampingParameter = 1;
    end
    obj.name = name;
    obj.numberOfPreSmoothingIterations = numberOfPreSmoothingIterations;
    obj.numberOfPostSmoothingIterations = numberOfPostSmoothingIterations;
    obj.dampingParameter = options.dampingParameter;
end

function approx = apply(obj,A,b,approx,j,type)
    arguments
        obj
        A
        b
        approx
        j % number current level
        type {mustBeMember(type,{'pre','post'})}
    end
    
    switch(type)
        case 'pre'
            if (length(obj.numberOfPreSmoothingIterations)>1)
                numberOfIterations = obj.numberOfPreSmoothingIterations(j);
            else
                numberOfIterations = obj.numberOfPreSmoothingIterations;
            end
        case 'post'
            if (length(obj.numberOfPostSmoothingIterations)>1)
                numberOfIterations = obj.numberOfPostSmoothingIterations(j);
            else
                numberOfIterations = obj.numberOfPostSmoothingIterations;
            end
    end
    
    switch(obj.name)
        case 'jacobi'
            approx = jacobi(A,b,approx,numberOfIterations,obj.dampingParameter);
        case 'gs'
            approx = gs(A,b,approx,numberOfIterations);
    end
end
end
end