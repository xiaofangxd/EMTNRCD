clc
clear %1,2,3,5,6,7,8,10
for i = 1:20
    main('-algorithm',@CommunitySparseEA,'-problem',@NREG1,'-evaluation',2e5,'-run',i)
    main('-algorithm',@CommunitySparseEA,'-problem',@NREG2,'-evaluation',2e5,'-run',i)
    main('-algorithm',@CommunitySparseEA,'-problem',@NREG3,'-evaluation',2e5,'-run',i)
    main('-algorithm',@CommunitySparseEA,'-problem',@NREG5,'-evaluation',2e5,'-run',i)
    main('-algorithm',@CommunitySparseEA,'-problem',@NREG6,'-evaluation',2e5,'-run',i)
    main('-algorithm',@CommunitySparseEA,'-problem',@NREG7,'-evaluation',2e5,'-run',i)
    main('-algorithm',@CommunitySparseEA,'-problem',@NREG8,'-evaluation',2e5,'-run',i)
    main('-algorithm',@CommunitySparseEA,'-problem',@NREG10,'-evaluation',2e5,'-run',i)
    main('-algorithm',@CommunitySparseEA,'-problem',@NRRN1,'-evaluation',2e5,'-run',i)
    main('-algorithm',@CommunitySparseEA,'-problem',@NRRN2,'-evaluation',2e5,'-run',i)
    main('-algorithm',@CommunitySparseEA,'-problem',@NRRN3,'-evaluation',2e5,'-run',i)
    main('-algorithm',@CommunitySparseEA,'-problem',@NRRN5,'-evaluation',2e5,'-run',i)
    main('-algorithm',@CommunitySparseEA,'-problem',@NRRN6,'-evaluation',2e5,'-run',i)
    main('-algorithm',@CommunitySparseEA,'-problem',@NRRN7,'-evaluation',2e5,'-run',i)
    main('-algorithm',@CommunitySparseEA,'-problem',@NRRN8,'-evaluation',2e5,'-run',i)
    main('-algorithm',@CommunitySparseEA,'-problem',@NRRN10,'-evaluation',2e5,'-run',i)
end
