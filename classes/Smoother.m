classdef Smoother
% SMOOTHER class for calling smoothing routines

properties
    numberOfPreSmoothingIterations  % vector - possible to choose different number on each level
    numberOfPostSmoothingIterations % vector - possible to choose different number on each level
end

methods 
    function smoother = Smoother(numberOfPreSmoothingIterations,numberOfPostSmoothingIterations)
    % SMOOTHER Construct an instance of class Smoother

    arguments
        numberOfPreSmoothingIterations
        numberOfPostSmoothingIterations
    end
    smoother.numberOfPreSmoothingIterations = numberOfPreSmoothingIterations;
    smoother.numberOfPostSmoothingIterations = numberOfPostSmoothingIterations;
    end
    
    function numberOfIterations = getNumberOfIterations(smoother,level,type)
        switch(type)
            case 'pre'
                if (length(smoother.numberOfPreSmoothingIterations)>1)
                    numberOfIterations = smoother.numberOfPreSmoothingIterations(level);
                else
                    numberOfIterations = smoother.numberOfPreSmoothingIterations;
                end
            case 'post'
                if (length(smoother.numberOfPostSmoothingIterations)>1)
                    numberOfIterations = smoother.numberOfPostSmoothingIterations(level);
                else
                    numberOfIterations = smoother.numberOfPostSmoothingIterations;
                end
        end
    end
end

methods (Abstract) 
    approx = apply(smoother)
end

end % class