MODULE PARALLEL_PRECON
  USE PRECISION;      USE global_variables;  USE MP_INTERFACE;
  USE maths;          USE gather_scatter;
  USE PARALLEL_SUPPLEMENTARY_MATHS;
  USE OMP_LIB;
  CONTAINS
!================================================
!
! Partial assembly Additive Schwarz
! solve (block-Jacobi) used to solve
! /precondition the global system.
!
!
!================================================

!---
! PARAFEMS classic Partitioner
!---
SUBROUTINE PF_PARTITIONER1(nn_pp1, nn_pp2, nn_pp, num, nn, npess, numpes)
  IMPLICIT NONE
  INTEGER,   INTENT(IN)   :: nn, npess, numpes
  INTEGER,   INTENT(INOUT):: nn_pp1, nn_pp2, nn_pp, num

  nn_pp2 = 0;
  nn_pp1 = 0;
  nn_pp  = 0;
  num    = 0;

  nn_pp2 = nn/npess
  num    = nn - nn_pp2*npess
  nn_pp1 = nn_pp2
  IF(num/=0) nn_pp1 = nn_pp1 + 1
  IF((numpes <= num).OR.(num == 0))THEN
    nn_pp = nn_pp1;
  ELSE
    nn_pp = nn_pp2;
  ENDIF
  RETURN
ENDSUBROUTINE PF_PARTITIONER1


!---
! Calculate the process that owns
! a particular node/element
!---
PURE FUNCTION FIND_NODE_PROC1(GnodeID,nn_pp1,nn_pp2,num,npess) RESULT(procID)
  IMPLICIT NONE
  INTEGER, INTENT(IN):: GnodeID, npess, num, nn_pp1, nn_pp2;
  INTEGER            :: PID1, A, B;
  INTEGER            :: procID

  A = GnodeID - num*nn_pp1
  IF(A == 0) procID = num;
  IF(A <  0)THEN
    PID1 = GnodeID/nn_pp1
    B = GnodeID - PID1*nn_pp1;
    IF(B == 0) procID = PID1
    IF(B /= 0) procID = PID1 + 1;
  ENDIF
  IF(A >  0)THEN
    PID1 = A/nn_pp2;
    B = A - PID1*nn_pp2
    IF(B == 0) procID = PID1 + num
    IF(B /= 0) procID = PID1 + num + 1;
  ENDIF
ENDFUNCTION FIND_NODE_PROC1


!---
! Calculates the local numberingID from global
! numberingID of a particular node/element
!---
PURE FUNCTION GLOBAL_TO_LOCAL1(G_nodeID,nn_pp1,nn_pp2,num,npess) RESULT(L_nodeID)
  IMPLICIT NONE
  INTEGER, INTENT(IN):: G_nodeID, num, nn_pp1, nn_pp2, npess;
  INTEGER            :: L_nodeID, procID, threshhold;

  threshhold = num*nn_pp1;
  procID = FIND_NODE_PROC1(G_nodeID,nn_pp1,nn_pp2,num,npess)
  procID = procID - 1;
  L_nodeID = G_nodeID - procID*nn_pp1
  IF(num /= 0)THEN
    IF(G_nodeID > threshhold)THEN
      L_nodeID = G_nodeID - num*nn_pp1 - (procID - num)*nn_pp2
    ENDIF
  ENDIF
ENDFUNCTION GLOBAL_TO_LOCAL1

!----
! Returns a matrix of absolute values
!----
SUBROUTINE ABS_MAT1(ABSAMat, Amat, N)
  IMPLICIT NONE
  INTEGER                 :: I, J;
  INTEGER,   INTENT(IN)   :: N;
  REAL(iwp), INTENT(IN)   :: Amat(N,N);
  REAL(iwp), INTENT(INOUT):: ABSAMat(N,N);
  REAL(iwp), PARAMETER    :: zero = 0._iwp;
  ABSAMat = zero;
  DO I=1,N
    DO J=1,N
      ABSAMat(I,J)=DABS(Amat(I,J))
    ENDDO
  ENDDO
  RETURN
ENDSUBROUTINE ABS_MAT1

!-------------------------------
! ILU(0) factorisation
!-------------------------------
SUBROUTINE ILU_factorisation(Lower, Upper, Amat, N, tol)
  USE precision;
  IMPLICIT NONE
  INTEGER                 :: I, J, K, L, M;
  INTEGER,   INTENT(IN)   :: N;
  REAL(iwp), INTENT(IN)   :: Amat(N,N), tol;
  REAL(iwp), INTENT(INOUT):: Lower(N,N), Upper(N,N);
  REAL(iwp), PARAMETER    :: one = 1._iwp, zero = 0._iwp;
  REAL(iwp)               :: AmatABS(N,N), LUMat(N,N);
  INTEGER                 :: DropRuleMat(N,N);

  !Drop rule ILU(0)
  DropRuleMat=0;
  CALL ABS_MAT1(AmatABS, Amat, N)
  DO I=1,N
    AmatABS(I,:) = AmatABS(I,:)/MAXVAL(AmatABS(I,:));
    DO J=1,N
      IF(AmatABS(I,J)<tol) DropRuleMat(I,J) = 0;
      IF(AmatABS(I,J)>tol) DropRuleMat(I,J) = 1;
    ENDDO
  ENDDO

  !ILU(0) factorisation
  LUMat = AMat
  DO I=1,N
    Upper(I,I)=one;
  ENDDO
  DO I=2,N
    DO K=1,(I-1); IF(DropRuleMat(I,K)==1)THEN
      LUMat(I,K) = LUMat(I,K)/LUMat(K,K)
      DO J=(K+1),N; IF(DropRuleMat(I,J)==1)THEN
        LUMat(I,J) = LUMat(I,J) - LUMat(I,K)*LUMat(K,J)
      ENDIF; ENDDO
    ENDIF; ENDDO
  ENDDO

  !Separate out the factors
  Lower = zero;
  Upper = zero;
  DO I=1,N
    Upper(I,I:N) = LUMat(I,I:N);
    Lower(I,1:I) = LUMat(I,1:I);
    Upper(I,I)   = one;
  ENDDO
  RETURN
END SUBROUTINE ILU_factorisation

!-------------------------------
! Calculate the number of
! neighboring processes
!-------------------------------
SUBROUTINE CalcNumNeighboringProcs(nodProcs_pp, nSendProcs, nRecvProcs, gg_pp &
                                 , nod, nel, nn, procID, nprocs)
  IMPLICIT NONE
  INTEGER               :: Iel, I, J;
  INTEGER, INTENT(IN)   :: nod, nel, nn, procID, nprocs
  INTEGER, INTENT(IN)   :: gg_pp(nod, nel);
  INTEGER, INTENT(INOUT):: nodProcs_pp(nod,nel); !procs where Nodes belong
  INTEGER, INTENT(INOUT):: nSendProcs;           !number of procs toSend 
  INTEGER, INTENT(INOUT):: nRecvProcs;           !number of procs toRecieve 
  INTEGER, ALLOCATABLE  :: procCommTable(:), procCommTable_tmp(:);            
  INTEGER               :: nn_pp1, nn_pp2, nn_pp, num;

  CALL PF_PARTITIONER1(nn_pp1, nn_pp2, nn_pp, num, nn, nprocs, procID)

  !
  ! Calculate the proc-
  ! number of each node
  !
  DO Iel = 1,nel
    DO I = 1,nod
      nodProcs_pp(I,Iel) = FIND_NODE_PROC1(gg_pp(I,Iel),nn_pp1,nn_pp2,num,nprocs)
    ENDDO
  ENDDO

  !
  ! Mark out the procs that the
  ! local process recieves from
  !
  ALLOCATE(procCommTable(nprocs), procCommTable_tmp(nprocs));
  procCommTable = 0;
  DO Iel = 1,nel
    DO I = 1,nod
      IF(procID /= nodProcs_pp(I,Iel))THEN
        procCommTable(nodProcs_pp(I,Iel)) = 1;
      ENDIF
    ENDDO
  ENDDO

  !
  ! Calculate the number of recieves
  ! from the procCommTable
  !
  nRecvProcs=0
  DO I = 1,nprocs
    IF(procCommTable(I) /= 0) nRecvProcs = nRecvProcs + 1;
  ENDDO


  !
  ! Calculate the number
  ! of procs to send to
  ! by using MPI_allreduce
  !
  procCommTable_tmp=0;
  CALL MPI_ALLREDUCE(procCommTable,procCommTable_tmp,nprocs,MPI_INTEGER,MPI_SUM &
                   , MPI_COMM_WORLD,ier);
  nSendProcs = procCommTable_tmp(procID);

  DEALLOCATE(procCommTable, procCommTable_tmp)
  RETURN
END SUBROUTINE CalcNumNeighboringProcs

!-------------------------------
! Calculate the number of non local components
!-------------------------------
SUBROUTINE NON_LOCAL_MAT_SUM(nodProcs_pp, nSendProcs, nRecvProcs, gg_pp &
                           , nod, nel, nn, procID, nprocs)
  IMPLICIT NONE
  INTEGER               :: Iel, I, J;
  INTEGER, INTENT(IN)   :: nod, nel, nn, procID, nprocs
  INTEGER, INTENT(IN)   :: gg_pp(nod, nel);
  INTEGER, INTENT(INOUT):: nodProcs_pp(nod,nel); !procs where Nodes belong
  INTEGER, INTENT(INOUT):: nSendProcs;           !number of procs toSend 
  INTEGER, INTENT(INOUT):: nRecvProcs;           !number of procs toRecieve 
  INTEGER, ALLOCATABLE  :: procCommTable(:), procCommTable_tmp(:);            
  INTEGER               :: nn_pp1, nn_pp2, nn_pp, num;
  INTEGER               :: procNodes(:,:)
  INTEGER               :: procSendSizes(:,:)
  INTEGER :: nsends, nrecvs

  CALL PF_PARTITIONER1(nn_pp1, nn_pp2, nn_pp, num, nn, nprocs, procID)

  procSends(nsends,nn_pp1)
  procs(nsends);
  procSendSizes(nsends)
  procID2


  !
  ! Find the size of the sends in each case
  !
  procSends = 0;
  DO Iel = 1,nel
    DO I = 1,nod
      procID2 = FIND_NODE_PROC1(gg_pp(I,Iel),nn_pp1,nn_pp2,num,nprocs) 
      L = GLOBAL_TO_LOCAL1(gg_pp(J,Iel),nn_pp1,nn_pp2,num,nprocs)
      PROC_SEARCH:DO J = 1,nsends
        K = J;
        IF(procID2 = procs(J)) EXIT PROC_SEARCH;
      ENDDO PROC_SEARCH
      procSends(K,L) = 1;
    ENDDO
  ENDDO

  procSendSizes = 0;
  DO I = 1,nSends
    DO J = 1,nn_pp1
      IF(procSends(I,J) /= 0) procSendSizes(I) = procSendSizes(I)  + 1;
    ENDDO
  ENDDO

  

  DEALLOCATE(procCommTable, procCommTable_tmp)
  RETURN
END SUBROUTINE NON_LOCAL_MAT_SUM

!-------------------------------
! partially assemble factorise
! matrix
!-------------------------------
SUBROUTINE BLOCK_JACOBI_ILU(LInv, UInv, stork_pp, gg_pp, MASK, nn &
                          , ntots, neqs_pp, nel, nodof, nod, procID, nprocs)
  IMPLICIT NONE
  INTEGER                 :: Iel, I, J, K, L, M, N, P, Q, procID1, procID2;
  INTEGER,   INTENT(IN)   :: nn, ntots, neqs_pp, nel, nodof, nod, procID, nprocs;
  INTEGER,   INTENT(IN)   :: gg_pp(nod,nel), MASK(nodof,nod)
  REAL(iwp), INTENT(IN)   :: stork_pp(ntots,ntots,nel);
  REAL(iwp), INTENT(INOUT):: LInv(neqs_pp,neqs_pp), UInv(neqs_pp,neqs_pp)
  REAL(iwp), PARAMETER    :: one = 1._iwp, zero = 0._iwp, tol=1.0E-15_iwp;
  REAL(iwp), ALLOCATABLE  :: lower(:,:), Upper(:,:), k_local(:,:), diag(:), diag_tmp(:,:);
  INTEGER                 :: nn_pp1, nn_pp2, nn_pp, num;

  INTEGER, ALLOCATABLE    :: nodProcs(:,:);
  INTEGER                 :: nSends, nRecvs

  ALLOCATE(k_local(neqs_pp,neqs_pp), diag(neqs_pp), diag_tmp(ntots,nel));
  ALLOCATE(lower(neqs_pp,neqs_pp), Upper(neqs_pp,neqs_pp));


  ALLOCATE(nodProcs(nod,nel))
  CALL PF_PARTITIONER1(nn_pp1, nn_pp2, nn_pp, num, nn, nprocs, procID)
  CALL CalcNumNeighboringProcs(nodProcs, nSends, nRecvs, gg_pp, nod, nel, nn, procID, nprocs)
  WRITE(*,*) procID, nSends, nRecvs;
  DEALLOCATE(nodProcs)

  k_local = zero;
  DO Iel = 1,nel
    DO I = 1,nod
      K = GLOBAL_TO_LOCAL1(gg_pp(I,Iel),nn_pp1,nn_pp2,num,nprocs)
      procID1 = FIND_NODE_PROC1(gg_pp(I,Iel),nn_pp1,nn_pp2,num,nprocs)
      IF(procID == procID1)THEN
        DO J = 1,nod
          L = GLOBAL_TO_LOCAL1(gg_pp(J,Iel),nn_pp1,nn_pp2,num,nprocs)
          procID2 = FIND_NODE_PROC1(gg_pp(J,Iel),nn_pp1,nn_pp2,num,nprocs) 
          IF(procID == procID2)THEN
            DO P = 1,nodof
              M = (P-1)*nn_pp + K;
              DO Q = 1,nodof
                N = (Q-1)*nn_pp + L;
                k_local(M,N) = k_local(M,N) + stork_pp(MASK(P,I),MASK(Q,J),Iel);
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ENDIF
    ENDDO
  ENDDO

  !Calculate the Diagonal block
  diag = zero;
  diag_tmp = zero;
  DO Iel=1,nel
    DO I = 1,ntots 
      diag_tmp(I,Iel) = diag_tmp(I,Iel) + stork_pp(I,I,Iel);
    ENDDO
  ENDDO
  CALL SCATTERM(diag, diag_tmp, MASK, ntots, nodof, nod, nel, neqs_pp, nn_pp) 
  DO I=1,neqs_pp
    diag(I) = DMAX1(DABS(diag(I) ), MAXVAL(DABS(k_local(I,:) ) ) )
    IF(diag(I) < tol) diag(I) = one;
    k_local(I,:) = k_local(I,:)/diag(I);
  ENDDO

  !Calculate the ILU Factorisation
  CALL ILU_factorisation(Lower, Upper, k_local, neqs_pp, tol)

  !Invert the matrices
  CALL Invert2(Lower, LInv, neqs_pp)
  CALL Invert2(Upper, UInv, neqs_pp)

  DO I=1,neqs_pp
    LInv(I,:) = LInv(I,:)/diag(I);
  ENDDO
  LInv = MATMUL(UInv,LInv)


  IF(numpe==1)THEN
    OPEN(76,FILE="Results2/ILUFactorsL.dat",STATUS='REPLACE',ACTION='WRITE')
    OPEN(77,FILE="Results2/ILUFactorsU.dat",STATUS='REPLACE',ACTION='WRITE')
    OPEN(78,FILE="Results2/ILUFactorsLInv.dat",STATUS='REPLACE',ACTION='WRITE')
    OPEN(79,FILE="Results2/ILUFactorsUInv.dat",STATUS='REPLACE',ACTION='WRITE')
    DO I=1,neqs_pp
      WRITE(76,*) Lower(I,:)
      WRITE(77,*) Upper(I,:)
      WRITE(78,*) LInv(I,:)
      WRITE(79,*) UInv(I,:)
    ENDDO
    CLOSE(76)
    CLOSE(77)
    CLOSE(78)
    CLOSE(79)
  ENDIF
  DEALLOCATE(Lower, Upper, k_local, diag, diag_tmp)
  RETURN
END SUBROUTINE BLOCK_JACOBI_ILU

!-------------------------------
! Matrix vector operator
!-------------------------------
SUBROUTINE MATVEC_Local(Ax_vec, Amat, x_vec, neqs_pp)
  INTEGER                 :: I, J;
  INTEGER                 :: nthreads, threadID, tnstart, tnend; !OpenMP stuff
  INTEGER,   INTENT(IN)   :: neqs_pp
  REAL(iwp), INTENT(IN)   :: Amat(neqs_pp,neqs_pp), x_vec(neqs_pp);
  REAL(iwp), INTENT(INOUT):: Ax_vec(neqs_pp); 

  Ax_vec = 0._iwp;
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I, J, threadID, nthreads, tnstart, tnend)
  !$OMP BARRIER
  threadID = OMP_GET_THREAD_NUM();
  nthreads = OMP_GET_MAX_THREADS();
  tnstart = ITERATOR_START(threadID, nthreads, neqs_pp)
  tnend   = ITERATOR_END(threadID, nthreads, neqs_pp) 

  DO I = tnstart,tnend
    DO J = 1,neqs_pp
      Ax_vec(I) = Ax_vec(I) + Amat(I,J)*x_vec(J)
    ENDDO
  ENDDO
  !$OMP BARRIER
  !$OMP END PARALLEL
  RETURN
END SUBROUTINE MATVEC_Local
!================================================
!================================================
!================================================
ENDMODULE PARALLEL_PRECON
