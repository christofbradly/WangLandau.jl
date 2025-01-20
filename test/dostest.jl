using BenchmarkTools

function logtest(a, b)
    a > b && return true
    return log(rand()) < a - b
end

function exptest(a, b)
    a > b && return true
    return rand() < exp(a) / exp(b)
end

function exptest2(a, b)
    a > b && return true
    return rand() < exp(a - b)
end

function exptest3(a, b)
    a > b && return true
    return @fastmath rand() < exp(a - b)
end

a = (1.0 + 0.5*rand()) * 1e11
b = (1.5 + 0.5rand()) * 1e11

@btime logtest($a, $b)
@btime exptest($a, $b)
@btime exptest2($a, $b)
@btime exptest3($a, $b)

a = (1.5 + 0.5rand()) * 1e11
b = (1.0 + 0.5rand()) * 1e11

@btime logtest($a, $b)
@btime exptest($a, $b)
@btime exptest2($a, $b)
@btime exptest3($a, $b)
