classdef GSSmoother < Smoother

    methods
        function gsSmoother = GSSmoother(numberOfPreSmoothingIterations,numberOfPostSmoothingIterations)

            arguments
                numberOfPreSmoothingIterations
                numberOfPostSmoothingIterations
            end
            gsSmoother@Smoother(numberOfPreSmoothingIterations,numberOfPostSmoothingIterations)
        end

        function approx = apply(gsSmoother,A,b,approx,level,type)
            arguments
                gsSmoother,A,b,approx,level
                type {mustBeMember(type,{'pre','post'})}
            end
            numberOfIterations = gsSmoother.getNumberOfIterations(level,type);
            
            approx = gs(A,b,approx,numberOfIterations);
        end
    end
end