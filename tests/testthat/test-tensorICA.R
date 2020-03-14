context("tensorICA")

test_that("# The size of Source tensor in output should not be larger than input tensor", {
    
    data(buccalbloodtensor)
    Dim.l <- EstDim(buccalbloodtensor$data)
    dim <- Dim.l$dim
    dJ <- Dim.l$dJ
    tica.o <- DoTICA(Data = buccalbloodtensor$data, dim = dim, method = "FOBI")
    dim1 <- dim(tica.o$TPCA.S)[1]
    dim2 <- dim(buccalbloodtensor$data)[1]
    expect_lte(dim1,dim2)
    dim1 <- dim(tica.o$TPCA.S)[2]
    dim2 <- dim(buccalbloodtensor$data)[2]
    expect_lte(dim1,dim2)
    dim1 <- dim(tica.o$TPCA.S)[3]
    dim2 <- dim(buccalbloodtensor$data)[3]
    expect_lte(dim1,dim2)
    
})

