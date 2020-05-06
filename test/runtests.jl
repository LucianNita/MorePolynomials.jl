using MorePolynomials
using Intervals
using Test



@testset "Global.jl" begin
    # generate psudorandom values for this
    checkinfsTestFunc(a,b,c,d) = checkinfs(Interval(c,d),Interval(a,b)) # a b local -> c d global
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

    global2localTestFunc(a,b,c,d,x) = global2local(Interval(a,b),Interval(c,d),x)# a b local -> c d global
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
