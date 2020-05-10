using MorePolynomials
using Intervals
using Test



@testset "Global.jl" begin
    x = LGRPoly([1.0,2,3])
    y = GlobalPoly(x,0,1)

    @test_throws ErrorException push!(y,x,1,2)
    @test_throws ErrorException push!(y,x,0,1)
    @test_throws ErrorException push!(y,x,-1,0)
    @test_nowarn push!(y,x,2,3)


    x = LGRPoly([1.0,2,1])
    y = GlobalPoly(x,0,1)
    @test_nowarn update!(y,[1.0,1,1.0])
    @test_nowarn push!(y,x,2,3)
    @test_nowarn push!(y,x,1,2)
    @test_nowarn push!(y,x,3,4)
    @test_nowarn push!(y,x,-1,0)
    @test_throws DomainError y(8) 
    x1 = LGRPoly([0.0,1,0])
    x2 = LGRPoly([0.0,2,0])
    x3 = LGRPoly([0.0,3,0])
    y1 = GlobalPoly(x1,0,1)
    push!(y1,x2,1,2)
    push!(y1,x3,2,3)
    @test_throws DimensionMismatch update!(y1,[1.0,1,1])
    @test_nowarn update!(y1,[0.0,1,2,3,4,5,6])

    x1 = LGRPoly([0.0,1,0])
    x2 = LGRPoly([0.0,2,0])
    x3 = LGRPoly([0.0,3,0])
    y1 = GlobalPoly(x1,0,1)
    push!(y1,x2,2,4)
    push!(y1,x3,4,5)
    @test_nowarn update!(y1,[0.0,1,2,4,5,6,7,8])

    # generate psudorandom values for this
    checkinfsTestFunc(a,b,c,d) = MorePolynomials.checkinfs(Interval(c,d),Interval(a,b)) # a b local -> c d global
    @test checkinfsTestFunc(-Inf,Inf,-Inf,Inf) == true
    @test_throws DomainError checkinfsTestFunc(-Inf,-2,-Inf,Inf) 
    @test checkinfsTestFunc(-2,2,-Inf,Inf) == true
    @test checkinfsTestFunc(-Inf,Inf,-Inf,2) == true
    @test checkinfsTestFunc(-Inf,Inf,-2,Inf) == true
    @test checkinfsTestFunc(-Inf,Inf,-2,2) == true
    @test checkinfsTestFunc(-Inf,2,-Inf,1) == true
    @test checkinfsTestFunc(-Inf,2,-Inf,3) == true
    @test_throws DomainError checkinfsTestFunc(-Inf,2,2,Inf)
    @test_throws DomainError checkinfsTestFunc(-2, Inf,-Inf,2)
    @test checkinfsTestFunc(2,Inf,3,Inf) == true
    @test checkinfsTestFunc(-2,2,3,Inf) == true
    @test checkinfsTestFunc(3,4,3,Inf) == true
    @test checkinfsTestFunc(3,4,-Inf,2) == true
    @test checkinfsTestFunc(3,Inf,3,4) == true
    @test_throws DomainError checkinfsTestFunc(3,Inf,2,4)
    @test_throws DomainError checkinfsTestFunc(3,Inf,1,2)
    @test checkinfsTestFunc(3,4,5,6) == true

    global2localTestFunc(a,b,c,d,x) = MorePolynomials.global2local(Interval(a,b),Interval(c,d),x)# a b local -> c d global
    @test global2localTestFunc(-Inf,Inf,-Inf,Inf,4) == 4
    @test global2localTestFunc(-2,2,-Inf,Inf,Inf) == 2
    @test global2localTestFunc(-2,2,-Inf,Inf,-Inf) == -2
    @test global2localTestFunc(-2,2,-Inf,Inf,0) == 0
    @test global2localTestFunc(-Inf,Inf,-Inf,2,1) == 1
    @test global2localTestFunc(-Inf,Inf,-2,Inf,1) == 1
    @test global2localTestFunc(-Inf,Inf,-2,2,1) == 1
    @test global2localTestFunc(2,Inf,3,Inf,3) == 2
    @test global2localTestFunc(2,Inf,3,Inf,Inf) == Inf
    @test global2localTestFunc(-Inf,2,-Inf,3,3) == 2
    @test global2localTestFunc(-Inf,2,-Inf,3,-Inf) == -Inf
    @test global2localTestFunc(3,4,-Inf,2,-Inf) == 3
    @test global2localTestFunc(3,4,2,Inf,2) == 3
    @test global2localTestFunc(3,4,2,Inf,Inf) == 4
    @test global2localTestFunc(3,4,-Inf,2,2) == 4
    @test global2localTestFunc(3,Inf,3,4,3) == 3
    @test global2localTestFunc(3,Inf,3,4,4) == 4
    @test global2localTestFunc(3,4,5,6,5) == 3
    @test global2localTestFunc(3,4,5,6,6) == 4
    @test global2localTestFunc(3,4,5,6,5.5) == 3.5
    

end
