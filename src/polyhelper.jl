# polyhelper function guesses what type of poly you want
export polyhelper


function polyhelper(y)
    LGRPoly(y)
end


function polyhelper(x,y)
    if length(x) < 6
        vandpoly(x,y)
    else
        LagrangePoly(x,y) # test speed of this against fitting function of Polynomials package
    end
end

#function polyhelper(x,y,b::Bound)
#    if length(x) < 6
#        vandpoly(x,y,b)
#    else
#        LagrangePoly(x,y,b)
#    end
#end


# add automatically gen GlobalPoly of local polys
