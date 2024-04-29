classdef JacobiSmoother < Smoother

    properties
        dampingParameter
    end
    methods
        function jacobiSmoother = JacobiSmoother(numberOfPreSmoothingIterations,numberOfPostSmoothingIterations,options)

            arguments
                numberOfPreSmoothingIterations
                numberOfPostSmoothingIterations
                options.dampingParameter = 1;
            end
            jacobiSmoother@Smoother(numberOfPreSmoothingIterations,numberOfPostSmoothingIterations)
            jacobiSmoother.dampingParameter = options.dampingParameter;
        end

        function approx = apply(jacobiSmoother,A,b,approx,level,type)
            arguments
                jacobiSmoother,A,b,approx,level
                type {mustBeMember(type,{'pre','post'})}
            end
            numberOfIterations = jacobiSmoother.getNumberOfIterations(level,type);
            
            approx = jacobi(A,b,approx,numberOfIterations,jacobiSmoother.dampingParameter);
        end
    end
end