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
unit rotations;
interface
uses Math, Ap, Sysutils;

procedure ApplyRotationsFromTheLeft(IsForward : Boolean;
     M1 : Integer;
     M2 : Integer;
     N1 : Integer;
     N2 : Integer;
     const C : TReal1DArray;
     const S : TReal1DArray;
     var A : TReal2DArray;
     var WORK : TReal1DArray);
procedure ApplyRotationsFromTheRight(IsForward : Boolean;
     M1 : Integer;
     M2 : Integer;
     N1 : Integer;
     N2 : Integer;
     const C : TReal1DArray;
     const S : TReal1DArray;
     var A : TReal2DArray;
     var WORK : TReal1DArray);
procedure GenerateRotation(F : Double;
     G : Double;
     var CS : Double;
     var SN : Double;
     var R : Double);

implementation

procedure TestRotations();forward;


(*************************************************************************
Application of a sequence of  elementary rotations to a matrix

The algorithm pre-multiplies the matrix by a sequence of rotation
transformations which is given by arrays C and S. Depending on the value
of the IsForward parameter either 1 and 2, 3 and 4 and so on (if IsForward=true)
rows are rotated, or the rows N and N-1, N-2 and N-3 and so on, are rotated.

Not the whole matrix but only a part of it is transformed (rows from M1 to
M2, columns from N1 to N2). Only the elements of this submatrix are changed.

Input parameters:
    IsForward   -   the sequence of the rotation application.
    M1,M2       -   the range of rows to be transformed.
    N1, N2      -   the range of columns to be transformed.
    C,S         -   transformation coefficients.
                    Array whose index ranges within [1..M2-M1].
    A           -   processed matrix.
    WORK        -   working array whose index ranges within [N1..N2].

Output parameters:
    A           -   transformed matrix.

Utility subroutine.
*************************************************************************)
procedure ApplyRotationsFromTheLeft(IsForward : Boolean;
     M1 : Integer;
     M2 : Integer;
     N1 : Integer;
     N2 : Integer;
     const C : TReal1DArray;
     const S : TReal1DArray;
     var A : TReal2DArray;
     var WORK : TReal1DArray);
var
    J : Integer;
    JP1 : Integer;
    CTEMP : Double;
    STEMP : Double;
    TEMP : Double;
begin
    if (M1>M2) or (N1>N2) then
    begin
        Exit;
    end;
    
    //
    // Form  P * A
    //
    if IsForward then
    begin
        if N1<>N2 then
        begin
            
            //
            // Common case: N1<>N2
            //
            J:=M1;
            while J<=M2-1 do
            begin
                CTEMP := C[J-M1+1];
                STEMP := S[J-M1+1];
                if (CTEMP<>1) or (STEMP<>0) then
                begin
                    JP1 := J+1;
                    APVMove(@WORK[0], N1, N2, @A[JP1][0], N1, N2, CTEMP);
                    APVSub(@WORK[0], N1, N2, @A[J][0], N1, N2, STEMP);
                    APVMul(@A[J][0], N1, N2, CTEMP);
                    APVAdd(@A[J][0], N1, N2, @A[JP1][0], N1, N2, STEMP);
                    APVMove(@A[JP1][0], N1, N2, @WORK[0], N1, N2);
                end;
                Inc(J);
            end;
        end
        else
        begin
            
            //
            // Special case: N1=N2
            //
            J:=M1;
            while J<=M2-1 do
            begin
                CTEMP := C[J-M1+1];
                STEMP := S[J-M1+1];
                if (CTEMP<>1) or (STEMP<>0) then
                begin
                    TEMP := A[J+1,N1];
                    A[J+1,N1] := CTEMP*TEMP-STEMP*A[J,N1];
                    A[J,N1] := STEMP*TEMP+CTEMP*A[J,N1];
                end;
                Inc(J);
            end;
        end;
    end
    else
    begin
        if N1<>N2 then
        begin
            
            //
            // Common case: N1<>N2
            //
            J:=M2-1;
            while J>=M1 do
            begin
                CTEMP := C[J-M1+1];
                STEMP := S[J-M1+1];
                if (CTEMP<>1) or (STEMP<>0) then
                begin
                    JP1 := J+1;
                    APVMove(@WORK[0], N1, N2, @A[JP1][0], N1, N2, CTEMP);
                    APVSub(@WORK[0], N1, N2, @A[J][0], N1, N2, STEMP);
                    APVMul(@A[J][0], N1, N2, CTEMP);
                    APVAdd(@A[J][0], N1, N2, @A[JP1][0], N1, N2, STEMP);
                    APVMove(@A[JP1][0], N1, N2, @WORK[0], N1, N2);
                end;
                Dec(J);
            end;
        end
        else
        begin
            
            //
            // Special case: N1=N2
            //
            J:=M2-1;
            while J>=M1 do
            begin
                CTEMP := C[J-M1+1];
                STEMP := S[J-M1+1];
                if (CTEMP<>1) or (STEMP<>0) then
                begin
                    TEMP := A[J+1,N1];
                    A[J+1,N1] := CTEMP*TEMP-STEMP*A[J,N1];
                    A[J,N1] := STEMP*TEMP+CTEMP*A[J,N1];
                end;
                Dec(J);
            end;
        end;
    end;
end;


(*************************************************************************
Application of a sequence of  elementary rotations to a matrix

The algorithm post-multiplies the matrix by a sequence of rotation
transformations which is given by arrays C and S. Depending on the value
of the IsForward parameter either 1 and 2, 3 and 4 and so on (if IsForward=true)
rows are rotated, or the rows N and N-1, N-2 and N-3 and so on are rotated.

Not the whole matrix but only a part of it is transformed (rows from M1
to M2, columns from N1 to N2). Only the elements of this submatrix are changed.

Input parameters:
    IsForward   -   the sequence of the rotation application.
    M1,M2       -   the range of rows to be transformed.
    N1, N2      -   the range of columns to be transformed.
    C,S         -   transformation coefficients.
                    Array whose index ranges within [1..N2-N1].
    A           -   processed matrix.
    WORK        -   working array whose index ranges within [M1..M2].

Output parameters:
    A           -   transformed matrix.

Utility subroutine.
*************************************************************************)
procedure ApplyRotationsFromTheRight(IsForward : Boolean;
     M1 : Integer;
     M2 : Integer;
     N1 : Integer;
     N2 : Integer;
     const C : TReal1DArray;
     const S : TReal1DArray;
     var A : TReal2DArray;
     var WORK : TReal1DArray);
var
    J : Integer;
    JP1 : Integer;
    CTEMP : Double;
    STEMP : Double;
    TEMP : Double;
    i_ : Integer;
begin
    
    //
    // Form A * P'
    //
    if IsForward then
    begin
        if M1<>M2 then
        begin
            
            //
            // Common case: M1<>M2
            //
            J:=N1;
            while J<=N2-1 do
            begin
                CTEMP := C[J-N1+1];
                STEMP := S[J-N1+1];
                if (CTEMP<>1) or (STEMP<>0) then
                begin
                    JP1 := J+1;
                    for i_ := M1 to M2 do
                    begin
                        WORK[i_] := CTEMP*A[i_,JP1];
                    end;
                    for i_ := M1 to M2 do
                    begin
                        WORK[i_] := WORK[i_] - STEMP*A[i_,J];
                    end;
                    for i_ := M1 to M2 do
                    begin
                        A[i_,J] := CTEMP*A[i_,J];
                    end;
                    for i_ := M1 to M2 do
                    begin
                        A[i_,J] := A[i_,J] + STEMP*A[i_,JP1];
                    end;
                    for i_ := M1 to M2 do
                    begin
                        A[i_,JP1] := WORK[i_];
                    end;
                end;
                Inc(J);
            end;
        end
        else
        begin
            
            //
            // Special case: M1=M2
            //
            J:=N1;
            while J<=N2-1 do
            begin
                CTEMP := C[J-N1+1];
                STEMP := S[J-N1+1];
                if (CTEMP<>1) or (STEMP<>0) then
                begin
                    TEMP := A[M1,J+1];
                    A[M1,J+1] := CTEMP*TEMP-STEMP*A[M1,J];
                    A[M1,J] := STEMP*TEMP+CTEMP*A[M1,J];
                end;
                Inc(J);
            end;
        end;
    end
    else
    begin
        if M1<>M2 then
        begin
            
            //
            // Common case: M1<>M2
            //
            J:=N2-1;
            while J>=N1 do
            begin
                CTEMP := C[J-N1+1];
                STEMP := S[J-N1+1];
                if (CTEMP<>1) or (STEMP<>0) then
                begin
                    JP1 := J+1;
                    for i_ := M1 to M2 do
                    begin
                        WORK[i_] := CTEMP*A[i_,JP1];
                    end;
                    for i_ := M1 to M2 do
                    begin
                        WORK[i_] := WORK[i_] - STEMP*A[i_,J];
                    end;
                    for i_ := M1 to M2 do
                    begin
                        A[i_,J] := CTEMP*A[i_,J];
                    end;
                    for i_ := M1 to M2 do
                    begin
                        A[i_,J] := A[i_,J] + STEMP*A[i_,JP1];
                    end;
                    for i_ := M1 to M2 do
                    begin
                        A[i_,JP1] := WORK[i_];
                    end;
                end;
                Dec(J);
            end;
        end
        else
        begin
            
            //
            // Special case: M1=M2
            //
            J:=N2-1;
            while J>=N1 do
            begin
                CTEMP := C[J-N1+1];
                STEMP := S[J-N1+1];
                if (CTEMP<>1) or (STEMP<>0) then
                begin
                    TEMP := A[M1,J+1];
                    A[M1,J+1] := CTEMP*TEMP-STEMP*A[M1,J];
                    A[M1,J] := STEMP*TEMP+CTEMP*A[M1,J];
                end;
                Dec(J);
            end;
        end;
    end;
end;


(*************************************************************************
The subroutine generates the elementary rotation, so that:

[  CS  SN  ]  .  [ F ]  =  [ R ]
[ -SN  CS  ]     [ G ]     [ 0 ]

CS**2 + SN**2 = 1
*************************************************************************)
procedure GenerateRotation(F : Double;
     G : Double;
     var CS : Double;
     var SN : Double;
     var R : Double);
var
    F1 : Double;
    G1 : Double;
begin
    if G=0 then
    begin
        CS := 1;
        SN := 0;
        R := F;
    end
    else
    begin
        if F=0 then
        begin
            CS := 0;
            SN := 1;
            R := G;
        end
        else
        begin
            F1 := F;
            G1 := G;
            R := SQRT(Sqr(F1)+Sqr(G1));
            CS := F1/R;
            SN := G1/R;
            if (ABSReal(F)>ABSReal(G)) and (CS<0) then
            begin
                CS := -CS;
                SN := -SN;
                R := -R;
            end;
        end;
    end;
end;


procedure TestRotations();
var
    AL1 : TReal2DArray;
    AL2 : TReal2DArray;
    AR1 : TReal2DArray;
    AR2 : TReal2DArray;
    CL : TReal1DArray;
    SL : TReal1DArray;
    CR : TReal1DArray;
    SR : TReal1DArray;
    W : TReal1DArray;
    M : Integer;
    N : Integer;
    MaxMN : Integer;
    T : Double;
    Pass : Integer;
    PassCount : Integer;
    I : Integer;
    J : Integer;
    Err : Double;
    MaxErr : Double;
    IsForward : Boolean;
begin
    PassCount := 1000;
    MaxErr := 0;
    Pass:=1;
    while Pass<=PassCount do
    begin
        
        //
        // settings
        //
        M := 2+RandomInteger(50);
        N := 2+RandomInteger(50);
        IsForward := RandomReal>0.5;
        MaxMN := Max(M, N);
        SetLength(AL1, M+1, N+1);
        SetLength(AL2, M+1, N+1);
        SetLength(AR1, M+1, N+1);
        SetLength(AR2, M+1, N+1);
        SetLength(CL, M-1+1);
        SetLength(SL, M-1+1);
        SetLength(CR, N-1+1);
        SetLength(SR, N-1+1);
        SetLength(W, MaxMN+1);
        
        //
        // matrices and rotaions
        //
        I:=1;
        while I<=M do
        begin
            J:=1;
            while J<=N do
            begin
                AL1[I,J] := 2*RandomReal-1;
                AL2[I,J] := AL1[I,J];
                AR1[I,J] := AL1[I,J];
                AR2[I,J] := AL1[I,J];
                Inc(J);
            end;
            Inc(I);
        end;
        I:=1;
        while I<=M-1 do
        begin
            T := 2*Pi*RandomReal;
            CL[I] := Cos(T);
            SL[I] := Sin(T);
            Inc(I);
        end;
        J:=1;
        while J<=N-1 do
        begin
            T := 2*Pi*RandomReal;
            CR[J] := Cos(T);
            SR[J] := Sin(T);
            Inc(J);
        end;
        
        //
        // Test left
        //
        ApplyRotationsFromTheLeft(IsForward, 1, M, 1, N, CL, SL, AL1, W);
        J:=1;
        while J<=N do
        begin
            ApplyRotationsFromTheLeft(IsForward, 1, M, J, J, CL, SL, AL2, W);
            Inc(J);
        end;
        Err := 0;
        I:=1;
        while I<=M do
        begin
            J:=1;
            while J<=N do
            begin
                Err := Max(Err, AbsReal(AL1[I,J]-AL2[I,J]));
                Inc(J);
            end;
            Inc(I);
        end;
        MaxErr := Max(Err, MaxErr);
        
        //
        // Test right
        //
        ApplyRotationsFromTheRight(IsForward, 1, M, 1, N, CR, SR, AR1, W);
        I:=1;
        while I<=M do
        begin
            ApplyRotationsFromTheRight(IsForward, I, I, 1, N, CR, SR, AR2, W);
            Inc(I);
        end;
        Err := 0;
        I:=1;
        while I<=M do
        begin
            J:=1;
            while J<=N do
            begin
                Err := Max(Err, AbsReal(AR1[I,J]-AR2[I,J]));
                Inc(J);
            end;
            Inc(I);
        end;
        MaxErr := Max(Err, MaxErr);
        Inc(Pass);
    end;
    Write(Format('TESTING ROTATIONS'#13#10'',[]));
    Write(Format('Pass count %0d'#13#10'',[
        PassCount]));
    Write(Format('Error is %5.4e'#13#10'',[
        MaxErr]));
end;


end.