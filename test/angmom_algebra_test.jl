@testset "GMP double factorial" begin
    for n in (-1, 0, 1)
        @test NuclearToolkit.mydoublefactorial(n) === 1.0
    end
    # Compare the native call with an independent exact product, including
    # results beyond machine-integer precision and Float64's finite range.
    for n in (2, 3, 4, 5, 20, 21, 64, 65, 170, 171, 200, 201, 300, 301, 302)
        expected = Float64(prod(big.(n:-2:1)))
        @test NuclearToolkit.mydoublefactorial(n) === expected
    end
    for n in (Int32(9), UInt32(10), big(21))
        @test NuclearToolkit.mydoublefactorial(n) === Float64(prod(big.(n:-2:1)))
    end
    @test_throws DomainError NuclearToolkit.mydoublefactorial(-2)
    @test_throws DomainError NuclearToolkit.mydoublefactorial(-3)
    # GMP takes unsigned long, which is narrower than UInt on 64-bit Windows.
    @test_throws InexactError NuclearToolkit.mydoublefactorial(big(typemax(Culong)) + 1)
end
