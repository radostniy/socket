(*************************************************************************
Copyright (c) 1992-2007 The University of Tennessee.  All rights reserved.

Contributors:
    * Sergey Bochkanov (ALGLIB project). Translation from FORTRAN to
      pseudocode.

See subroutines comments for additional copyrights.

>>> SOURCE LICENSE >>>
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation (www.fsf.org); either version 2 of the 
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License is available at
http://www.fsf.org/licensing/licenses

>>> END OF LICENSE >>>
*************************************************************************)
unit reflections;
interface
uses Math, Ap, Sysutils;

procedure GenerateReflection(var X : TReal1DArray;
     N : Integer;
     var Tau : Double);
procedure ApplyReflectionFromTheLeft(var C : TReal2DArray;
     Tau : Double;
     const V : TReal1DArray;
     M1 : Integer;
     M2 : Integer;
     N1 : Integer;
     N2 : Integer;
     var WORK : TReal1DArray);
procedure ApplyReflectionFromTheRight(var C : TReal2DArray;
     Tau : Double;
     const V : TReal1DArray;
     M1 : Integer;
     M2 : Integer;
     N1 : Integer;
     N2 : Integer;
     var WORK : TReal1DArray);

implementation

procedure TestReflections();forward;


(*************************************************************************
Generation of an elementary reflection transformation

The subroutine generates elementary reflection H of order N, so that, for
a given X, the following equality holds true:

    ( X(1) )   ( Beta )
H * (  ..  ) = (  0   )
    ( X(n) )   (  0   )

where
              ( V(1) )
H = 1 - Tau * (  ..  ) * ( V(1), ..., V(n) )
              ( V(n) )

where the first component of vector V equals 1.

Input parameters:
    X   -   vector. Array whose index ranges within [1..N].
    N   -   reflection order.

Output parameters:
    X   -   components from 2 to N are replaced with vector V.
            The first component is replaced with parameter Beta.
    Tau -   scalar value Tau. If X is a null vector, Tau equals 0,
            otherwise 1 <= Tau <= 2.

This subroutine is the modification of the DLARFG subroutines from
the LAPACK library. It has a similar functionality except for the
fact that it doesn’t handle errors when the intermediate results
cause an overflow.


MODIFICATIONS:
    24.12.2005 sign(Alpha) was replaced with an analogous to the Fortran SIGN code.

  -- LAPACK auxiliary routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     September 30, 1994
*************************************************************************)
procedure GenerateReflection(var X : TReal1DArray;
     N : Integer;
     var Tau : Double);
var
    J : Integer;
    Alpha : Double;
    XNORM : Double;
    V : Double;
    Beta : Double;
    MX : Double;
begin
    
    //
    // Executable Statements ..
    //
    if N<=1 then
    begin
        Tau := 0;
        Exit;
    end;
    
    //
    // XNORM = DNRM2( N-1, X, INCX )
    //
    Alpha := X[1];
    MX := 0;
    J:=2;
    while J<=N do
    begin
        MX := Max(AbsReal(X[J]), MX);
        Inc(J);
    end;
    XNORM := 0;
    if MX<>0 then
    begin
        J:=2;
        while J<=N do
        begin
            XNORM := XNORM+Sqr(X[J]/MX);
            Inc(J);
        end;
        XNORM := Sqrt(XNORM)*MX;
    end;
    if XNORM=0 then
    begin
        
        //
        // H  =  I
        //
        TAU := 0;
        Exit;
    end;
    
    //
    // general case
    //
    MX := Max(AbsReal(Alpha), AbsReal(XNORM));
    Beta := -MX*Sqrt(Sqr(Alpha/MX)+Sqr(XNORM/MX));
    if Alpha<0 then
    begin
        Beta := -Beta;
    end;
    TAU := (BETA-ALPHA)/BETA;
    V := 1/(Alpha-Beta);
    APVMul(@X[0], 2, N, V);
    X[1] := Beta;
end;


(*************************************************************************
Application of an elementary reflection to a rectangular matrix of size MxN

The algorithm pre-multiplies the matrix by an elementary reflection transformation
which is given by column V and scalar Tau (see the description of the
GenerateReflection procedure). Not the whole matrix but only a part of it
is transformed (rows from M1 to M2, columns from N1 to N2). Only the elements
of this submatrix are changed.

Input parameters:
    C       -   matrix to be transformed.
    Tau     -   scalar defining the transformation.
    V       -   column defining the transformation.
                Array whose index ranges within [1..M2-M1+1].
    M1, M2  -   range of rows to be transformed.
    N1, N2  -   range of columns to be transformed.
    WORK    -   working array whose indexes goes from N1 to N2.

Output parameters:
    C       -   the result of multiplying the input matrix C by the
                transformation matrix which is given by Tau and V.
                If N1>N2 or M1>M2, C is not modified.

  -- LAPACK auxiliary routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     September 30, 1994
*************************************************************************)
procedure ApplyReflectionFromTheLeft(var C : TReal2DArray;
     Tau : Double;
     const V : TReal1DArray;
     M1 : Integer;
     M2 : Integer;
     N1 : Integer;
     N2 : Integer;
     var WORK : TReal1DArray);
var
    T : Double;
    I : Integer;
    //VM : Integer;
begin
    if (Tau=0) or (N1>N2) or (M1>M2) then
    begin
        Exit;
    end;
    
    //
    // w := C' * v
    //
    //VM := M2-M1+1;
    I:=N1;
    while I<=N2 do
    begin
        WORK[I] := 0;
        Inc(I);
    end;
    I:=M1;
    while I<=M2 do
    begin
        T := V[I+1-M1];
        APVAdd(@WORK[0], N1, N2, @C[I][0], N1, N2, T);
        Inc(I);
    end;
    
    //
    // C := C - tau * v * w'
    //
    I:=M1;
    while I<=M2 do
    begin
        T := V[I-M1+1]*TAU;
        APVSub(@C[I][0], N1, N2, @WORK[0], N1, N2, T);
        Inc(I);
    end;
end;


(*************************************************************************
Application of an elementary reflection to a rectangular matrix of size MxN

The algorithm post-multiplies the matrix by an elementary reflection transformation
which is given by column V and scalar Tau (see the description of the
GenerateReflection procedure). Not the whole matrix but only a part of it
is transformed (rows from M1 to M2, columns from N1 to N2). Only the
elements of this submatrix are changed.

Input parameters:
    C       -   matrix to be transformed.
    Tau     -   scalar defining the transformation.
    V       -   column defining the transformation.
                Array whose index ranges within [1..N2-N1+1].
    M1, M2  -   range of rows to be transformed.
    N1, N2  -   range of columns to be transformed.
    WORK    -   working array whose indexes goes from M1 to M2.

Output parameters:
    C       -   the result of multiplying the input matrix C by the
                transformation matrix which is given by Tau and V.
                If N1>N2 or M1>M2, C is not modified.

  -- LAPACK auxiliary routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     September 30, 1994
*************************************************************************)
procedure ApplyReflectionFromTheRight(var C : TReal2DArray;
     Tau : Double;
     const V : TReal1DArray;
     M1 : Integer;
     M2 : Integer;
     N1 : Integer;
     N2 : Integer;
     var WORK : TReal1DArray);
var
    T : Double;
    I : Integer;
    VM : Integer;
begin
    if (Tau=0) or (N1>N2) or (M1>M2) then
    begin
        Exit;
    end;
    
    //
    // w := C * v
    //
    VM := N2-N1+1;
    I:=M1;
    while I<=M2 do
    begin
        T := APVDotProduct(@C[I][0], N1, N2, @V[0], 1, VM);
        WORK[I] := T;
        Inc(I);
    end;
    
    //
    // C := C - w * v'
    //
    I:=M1;
    while I<=M2 do
    begin
        T := WORK[I]*TAU;
        APVSub(@C[I][0], N1, N2, @V[0], 1, VM, T);
        Inc(I);
    end;
end;


procedure TestReflections();
var
    I : Integer;
    J : Integer;
    N : Integer;
    M : Integer;
    MaxMN : Integer;
    X : TReal1DArray;
    V : TReal1DArray;
    WORK : TReal1DArray;
    H : TReal2DArray;
    A : TReal2DArray;
    B : TReal2DArray;
    C : TReal2DArray;
    Tmp : Double;
    Beta : Double;
    Tau : Double;
    Err : Double;
    MER : Double;
    MEL : Double;
    MEG : Double;
    Pass : Integer;
    PassCount : Integer;
    i_ : Integer;
begin
    PassCount := 1000;
    MER := 0;
    MEL := 0;
    MEG := 0;
    Pass:=1;
    while Pass<=PassCount do
    begin
        
        //
        // Task
        //
        N := 1+RandomInteger(10);
        M := 1+RandomInteger(10);
        MaxMN := Max(M, N);
        
        //
        // Initialize
        //
        SetLength(X, MaxMN+1);
        SetLength(V, MaxMN+1);
        SetLength(WORK, MaxMN+1);
        SetLength(H, MaxMN+1, MaxMN+1);
        SetLength(A, MaxMN+1, MaxMN+1);
        SetLength(B, MaxMN+1, MaxMN+1);
        SetLength(C, MaxMN+1, MaxMN+1);
        
        //
        // GenerateReflection
        //
        I:=1;
        while I<=N do
        begin
            X[I] := 2*RandomReal-1;
            V[I] := X[I];
            Inc(I);
        end;
        GenerateReflection(V, N, Tau);
        Beta := V[1];
        V[1] := 1;
        I:=1;
        while I<=N do
        begin
            J:=1;
            while J<=N do
            begin
                if I=J then
                begin
                    H[I,J] := 1-Tau*V[I]*V[J];
                end
                else
                begin
                    H[I,J] := -Tau*V[I]*V[J];
                end;
                Inc(J);
            end;
            Inc(I);
        end;
        Err := 0;
        I:=1;
        while I<=N do
        begin
            Tmp := APVDotProduct(@H[I][0], 1, N, @X[0], 1, N);
            if I=1 then
            begin
                Err := Max(Err, AbsReal(Tmp-Beta));
            end
            else
            begin
                Err := Max(Err, AbsReal(Tmp));
            end;
            Inc(I);
        end;
        MEG := Max(MEG, Err);
        
        //
        // ApplyReflectionFromTheLeft
        //
        I:=1;
        while I<=M do
        begin
            X[I] := 2*RandomReal-1;
            V[I] := X[I];
            Inc(I);
        end;
        I:=1;
        while I<=M do
        begin
            J:=1;
            while J<=N do
            begin
                A[I,J] := 2*RandomReal-1;
                B[I,J] := A[I,J];
                Inc(J);
            end;
            Inc(I);
        end;
        GenerateReflection(V, M, Tau);
        //Beta := V[1];
        V[1] := 1;
        ApplyReflectionFromTheLeft(B, Tau, V, 1, M, 1, N, WORK);
        I:=1;
        while I<=M do
        begin
            J:=1;
            while J<=M do
            begin
                if I=J then
                begin
                    H[I,J] := 1-Tau*V[I]*V[J];
                end
                else
                begin
                    H[I,J] := -Tau*V[I]*V[J];
                end;
                Inc(J);
            end;
            Inc(I);
        end;
        I:=1;
        while I<=M do
        begin
            J:=1;
            while J<=N do
            begin
                Tmp := 0.0;
                for i_ := 1 to M do
                begin
                    Tmp := Tmp + H[I,i_]*A[i_,J];
                end;
                C[I,J] := Tmp;
                Inc(J);
            end;
            Inc(I);
        end;
        Err := 0;
        I:=1;
        while I<=M do
        begin
            J:=1;
            while J<=N do
            begin
                Err := Max(Err, AbsReal(B[I,J]-C[I,J]));
                Inc(J);
            end;
            Inc(I);
        end;
        MEL := Max(MEL, Err);
        
        //
        // ApplyReflectionFromTheRight
        //
        I:=1;
        while I<=N do
        begin
            X[I] := 2*RandomReal-1;
            V[I] := X[I];
            Inc(I);
        end;
        I:=1;
        while I<=M do
        begin
            J:=1;
            while J<=N do
            begin
                A[I,J] := 2*RandomReal-1;
                B[I,J] := A[I,J];
                Inc(J);
            end;
            Inc(I);
        end;
        GenerateReflection(V, N, Tau);
        //Beta := V[1];
        V[1] := 1;
        ApplyReflectionFromTheRight(B, Tau, V, 1, M, 1, N, WORK);
        I:=1;
        while I<=N do
        begin
            J:=1;
            while J<=N do
            begin
                if I=J then
                begin
                    H[I,J] := 1-Tau*V[I]*V[J];
                end
                else
                begin
                    H[I,J] := -Tau*V[I]*V[J];
                end;
                Inc(J);
            end;
            Inc(I);
        end;
        I:=1;
        while I<=M do
        begin
            J:=1;
            while J<=N do
            begin
                Tmp := 0.0;
                for i_ := 1 to N do
                begin
                    Tmp := Tmp + A[I,i_]*H[i_,J];
                end;
                C[I,J] := Tmp;
                Inc(J);
            end;
            Inc(I);
        end;
        Err := 0;
        I:=1;
        while I<=M do
        begin
            J:=1;
            while J<=N do
            begin
                Err := Max(Err, AbsReal(B[I,J]-C[I,J]));
                Inc(J);
            end;
            Inc(I);
        end;
        MER := Max(MER, Err);
        Inc(Pass);
    end;
    
    //
    // Overflow crash test
    //
    SetLength(X, 10+1);
    SetLength(V, 10+1);
    I:=1;
    while I<=10 do
    begin
        V[I] := MaxRealNumber*0.01*(2*RandomReal-1);
        Inc(I);
    end;
    GenerateReflection(V, 10, Tau);
    Write(Format('TESTING REFLECTIONS'#13#10'',[]));
    Write(Format('Pass count is %0d'#13#10'',[
        PassCount]));
    Write(Format('Generate     absolute error is       %5.4e'#13#10'',[
        MEG]));
    Write(Format('Apply(Left)  absolute error is       %5.4e'#13#10'',[
        MEL]));
    Write(Format('Apply(Right) absolute error is       %5.4e'#13#10'',[
        MER]));
    Write(Format('Overflow crash test passed'#13#10'',[]));
end;


end.