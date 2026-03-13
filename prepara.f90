! prepara.f90
!
! Reads sparg output (cadeias.txt, monomeros.txt) and writes the two
! input files required by conducao (vectorcad.txt, tabcad.txt).
!
! Parameters (dimx, dimy, dimz, inter, cross, aproxima) are read from
! conducao.ini using the same keyword format as conducao itself.
! If conducao.ini is absent, built-in defaults are used.
!
! All files must be in the working directory.
!
! Author: Ricardo Mendes Ribeiro

PROGRAM prepara

  IMPLICIT NONE

  ! Maximum chain length written by sparg (tam_max_cad in sparg/comum.f90).
  INTEGER, PARAMETER :: tam_max_cad = 20

  ! Maximum number of inter-chain connections per (chain, endpoint).
  ! Must be >= numeroligcad in conducao/comum module (currently 300).
  INTEGER, PARAMETER :: maxconn = 300

  ! ---- run parameters (read from conducao.ini) ----
  REAL(KIND=4) :: dimx_val, dimy_val, dimz_val
  REAL(KIND=4) :: inter, cross, aproxima, cutoff, cutoff2

  ! ---- chain data (from cadeias.txt) ----
  INTEGER :: Ncad
  INTEGER, ALLOCATABLE :: first_mon(:)   ! index of first monomer of each chain
  INTEGER, ALLOCATABLE :: ncad_chain(:)  ! number of monomeros per chain

  ! ---- monomer data (from monomeros.txt) ----
  INTEGER :: Nmon
  REAL(KIND=4), ALLOCATABLE :: px(:), py(:), pz(:)

  ! ---- chain geometry ----
  ! Vcad(i, 1/2/3) = start / middle / end position of chain i
  REAL(KIND=4), ALLOCATABLE :: Vcadx(:,:), Vcady(:,:), Vcadz(:,:)
  REAL(KIND=4), ALLOCATABLE :: Ccadx(:), Ccady(:), Ccadz(:) ! direction cosines

  ! ---- electrode chain counts ----
  INTEGER :: nbase, ntopo

  ! ---- inter-chain connectivity ----
  INTEGER, ALLOCATABLE :: nconn(:,:)             ! (Ncad,3) connections per endpoint
  INTEGER, ALLOCATABLE :: conn_dest(:,:,:)        ! (Ncad,3,maxconn) destination chain
  REAL(KIND=4), ALLOCATABLE :: conn_vec(:,:,:,:)  ! (Ncad,3,maxconn,4) dx,dy,dz,dist

  ! ---- work variables ----
  INTEGER  :: i, j, kk, mm, nc, mid_idx
  INTEGER  :: tmp_arr(tam_max_cad), idx_dum
  REAL(KIND=4) :: alf_dum, bet_dum, gam_dum
  REAL(KIND=4) :: x1, y1, z1, dx, dy, dz, d2
  REAL(KIND=4) :: vx, vy, vz, vlen
  REAL(KIND=4) :: best_d2, best_dx, best_dy, best_dz
  LOGICAL :: fexists

  ! ================================================================
  ! 1. Parameters
  ! ================================================================

  dimx_val = 20.0;  dimy_val = 20.0;  dimz_val = 100.0
  inter    =  0.65;  cross   =  0.15;  aproxima =   0.1

  INQUIRE(FILE='conducao.ini', EXIST=fexists)
  IF (fexists) THEN
    CALL read_ini(dimx_val, dimy_val, dimz_val, inter, cross, aproxima)
    WRITE(*,'(a)') ' Parameters from conducao.ini:'
  ELSE
    WRITE(*,'(a)') ' Warning: conducao.ini not found, using defaults.'
  ENDIF
  WRITE(*,'(a,3f8.2)') '   dimx, dimy, dimz (nm) = ', dimx_val, dimy_val, dimz_val
  WRITE(*,'(a,3f8.3)') '   inter, cross, aproxima (nm) = ', inter, cross, aproxima

  cutoff  = 2.5 * inter    ! include chains up to 2.5*inter; prob decays as exp(-(d-inter)/inter)
  cutoff2 = cutoff * cutoff

  ! ================================================================
  ! 2. Read cadeias.txt
  !
  ! Format (sparg print.f90):
  !   line 1:        Ncad
  !   lines 2..Ncad+1:  tam_max_cad integers; first is the index of the
  !                     first monomer, next ncad_chain(i)-1 are consecutive,
  !                     trailing zeros fill to tam_max_cad
  !   lines Ncad+2..:   ncadeia(i), one per line
  ! ================================================================

  INQUIRE(FILE='cadeias.txt', EXIST=fexists)
  IF (.NOT. fexists) THEN
    WRITE(*,*) 'Error: cadeias.txt not found.'
    STOP
  ENDIF

  OPEN(UNIT=1, FILE='cadeias.txt', STATUS='OLD')
  READ(1,*) Ncad
  ALLOCATE(first_mon(Ncad), ncad_chain(Ncad))
  DO i = 1, Ncad
    READ(1,*) tmp_arr            ! reads exactly tam_max_cad integers
    first_mon(i) = tmp_arr(1)   ! first monomer index of this chain
  ENDDO
  DO i = 1, Ncad
    READ(1,*) ncad_chain(i)
  ENDDO
  CLOSE(1)
  WRITE(*,'(a,i8)') ' Chains read: ', Ncad

  ! ================================================================
  ! 3. Read monomeros.txt
  !
  ! Format (sparg print.f90):
  !   line 1:  Nmon
  !   lines 2..: i  x  y  z  alfa  beta  gama    format (i7,6f12.5)
  ! ================================================================

  INQUIRE(FILE='monomeros.txt', EXIST=fexists)
  IF (.NOT. fexists) THEN
    WRITE(*,*) 'Error: monomeros.txt not found.'
    STOP
  ENDIF

  OPEN(UNIT=2, FILE='monomeros.txt', STATUS='OLD')
  READ(2,*) Nmon
  ALLOCATE(px(Nmon), py(Nmon), pz(Nmon))
  DO i = 1, Nmon
    READ(2,*) idx_dum, px(i), py(i), pz(i), alf_dum, bet_dum, gam_dum
  ENDDO
  CLOSE(2)
  WRITE(*,'(a,i8)') ' Monomers read: ', Nmon

  ! ================================================================
  ! 4. Compute chain geometry
  !
  ! Position 1 = start (first monomer)
  ! Position 2 = middle (monomer at index first_mon + (n-1)/2)
  ! Position 3 = end   (last monomer, first_mon + ncad_chain - 1)
  !
  ! Direction cosines: unit vector from start to end.
  ! ================================================================

  ALLOCATE(Vcadx(Ncad,3), Vcady(Ncad,3), Vcadz(Ncad,3))
  ALLOCATE(Ccadx(Ncad), Ccady(Ncad), Ccadz(Ncad))

  DO i = 1, Ncad
    ! start
    j = first_mon(i)
    Vcadx(i,1) = px(j);  Vcady(i,1) = py(j);  Vcadz(i,1) = pz(j)

    ! end
    j = first_mon(i) + ncad_chain(i) - 1
    Vcadx(i,3) = px(j);  Vcady(i,3) = py(j);  Vcadz(i,3) = pz(j)

    ! middle
    mid_idx = first_mon(i) + (ncad_chain(i) - 1) / 2
    Vcadx(i,2) = px(mid_idx)
    Vcady(i,2) = py(mid_idx)
    Vcadz(i,2) = pz(mid_idx)

    ! direction cosines: start -> end unit vector
    vx = Vcadx(i,3) - Vcadx(i,1)
    vy = Vcady(i,3) - Vcady(i,1)
    vz = Vcadz(i,3) - Vcadz(i,1)
    vlen = SQRT(vx*vx + vy*vy + vz*vz)
    IF (vlen > 0.) THEN
      Ccadx(i) = vx/vlen;  Ccady(i) = vy/vlen;  Ccadz(i) = vz/vlen
    ELSE
      Ccadx(i) = 0.;  Ccady(i) = 0.;  Ccadz(i) = 1.  ! fallback: along z
    ENDIF
  ENDDO

  ! ================================================================
  ! 5. Count base / topo chains
  !    (same criterion as conducao INPUTS subroutine)
  ! ================================================================

  nbase = 0;  ntopo = 0
  DO i = 1, Ncad
    IF      (Vcadz(i,1) <= aproxima .OR. Vcadz(i,3) <= aproxima) THEN
      nbase = nbase + 1
    ELSE IF (Vcadz(i,1) >= dimz_val - aproxima .OR. &
             Vcadz(i,3) >= dimz_val - aproxima) THEN
      ntopo = ntopo + 1
    ENDIF
  ENDDO
  WRITE(*,'(a,i6,a,i6)') ' Base chains: ', nbase, '   Top chains: ', ntopo

  ! ================================================================
  ! 6. Write vectorcad.txt
  !
  ! Expected by conducao SETPAR:
  !   line 1:        Ncad  base  topo          (free format)
  !   lines 2..Ncad+1:  ncadeia(i)  x1 y1 z1  x2 y2 z2  x3 y3 z3  cx cy cz
  !                                             format (i7,12f10.4)
  ! ================================================================

  OPEN(UNIT=3, FILE='vectorcad.txt')
  WRITE(3,*) Ncad, nbase, ntopo
  DO i = 1, Ncad
    WRITE(3,'(i7,12f10.4)') ncad_chain(i), &
      Vcadx(i,1), Vcady(i,1), Vcadz(i,1), &
      Vcadx(i,2), Vcady(i,2), Vcadz(i,2), &
      Vcadx(i,3), Vcady(i,3), Vcadz(i,3), &
      Ccadx(i),   Ccady(i),   Ccadz(i)
  ENDDO
  CLOSE(3)
  WRITE(*,*) 'vectorcad.txt written.'

  ! ================================================================
  ! 7. Compute inter-chain connectivity
  !
  ! For each source (chain i, endpoint kk=1,2,3) find all chains j
  ! within cutoff distance.  The jump vector points from the source
  ! position to the NEAREST of the three endpoints of chain j.
  !
  ! Minimum-image convention is applied in x and y (periodic BC
  ! matching conducao's POSICOES subroutine); z is open.
  ! ================================================================

  ALLOCATE(nconn(Ncad,3))
  ALLOCATE(conn_dest(Ncad,3,maxconn))
  ALLOCATE(conn_vec(Ncad,3,maxconn,4))
  nconn     = 0
  conn_dest = 0
  conn_vec  = 0.

  WRITE(*,'(a,f6.3,a)') ' Computing connectivity (cutoff =', cutoff, ' nm)...'

  DO i = 1, Ncad
    IF (MOD(i, MAX(1, Ncad/4)) == 0) &
      WRITE(*,'(a,i3,a)') '  ', 100*i/Ncad, '%'

    DO kk = 1, 3
      x1 = Vcadx(i,kk);  y1 = Vcady(i,kk);  z1 = Vcadz(i,kk)

      DO j = 1, Ncad
        IF (j == i) CYCLE

        ! nearest endpoint of chain j, with minimum-image in x,y
        best_d2 = 1.e30
        best_dx = 0.;  best_dy = 0.;  best_dz = 0.
        DO mm = 1, 3
          dx = Vcadx(j,mm) - x1
          dy = Vcady(j,mm) - y1
          dz = Vcadz(j,mm) - z1
          dx = dx - ANINT(dx/dimx_val)*dimx_val   ! nearest image x
          dy = dy - ANINT(dy/dimy_val)*dimy_val   ! nearest image y
          d2 = dx*dx + dy*dy + dz*dz
          IF (d2 < best_d2) THEN
            best_d2 = d2;  best_dx = dx;  best_dy = dy;  best_dz = dz
          ENDIF
        ENDDO

        IF (best_d2 > cutoff2) CYCLE

        nc = nconn(i,kk) + 1
        IF (nc > maxconn) THEN
          WRITE(*,'(a,2i6)') ' Warning: maxconn exceeded at chain,endpoint:', i, kk
          CYCLE
        ENDIF
        nconn(i,kk)         = nc
        conn_dest(i,kk,nc)  = j
        conn_vec(i,kk,nc,1) = best_dx
        conn_vec(i,kk,nc,2) = best_dy
        conn_vec(i,kk,nc,3) = best_dz
        conn_vec(i,kk,nc,4) = SQRT(best_d2)
      ENDDO
    ENDDO
  ENDDO
  WRITE(*,*) '  100%'

  ! ================================================================
  ! 8. Write tabcad.txt
  !
  ! Expected by conducao SETPAR:
  !   header:  ii  kk  jj             format (3i7)
  !            (chain, endpoint 1/2/3, number of jumps)
  !   per jump: dest  dx  dy  dz  dist  format (i7,4f10.4)
  !   (only blocks with jj > 0 are written)
  ! ================================================================

  OPEN(UNIT=4, FILE='tabcad.txt')
  DO i = 1, Ncad
    DO kk = 1, 3
      IF (nconn(i,kk) == 0) CYCLE
      WRITE(4,'(3i7)') i, kk, nconn(i,kk)
      DO nc = 1, nconn(i,kk)
        WRITE(4,'(i7,4f10.4)') conn_dest(i,kk,nc), &
          conn_vec(i,kk,nc,1), conn_vec(i,kk,nc,2), &
          conn_vec(i,kk,nc,3), conn_vec(i,kk,nc,4)
      ENDDO
    ENDDO
  ENDDO
  CLOSE(4)
  WRITE(*,*) 'tabcad.txt written.'

  STOP

CONTAINS

  !--------------------------------------------------------------------
  SUBROUTINE read_ini(dimx, dimy, dimz, inter, cross, aproxima)
  !--------------------------------------------------------------------
  ! Reads conducao.ini using the same keyword/value format as conducao.
  ! Only the four parameters needed here are extracted; everything else
  ! is skipped.
  !--------------------------------------------------------------------
    REAL(KIND=4), INTENT(INOUT) :: dimx, dimy, dimz, inter, cross, aproxima
    CHARACTER(25) :: wdkey
    CHARACTER(3)  :: key3
    INTEGER       :: ios_loc

    OPEN(UNIT=10, FILE='conducao.ini', STATUS='OLD')
    DO
      READ(10, '(A25)', IOSTAT=ios_loc) wdkey
      IF (ios_loc /= 0) EXIT
      key3 = wdkey(1:3)
      IF      (key3 == 'dim' .OR. key3 == 'DIM') THEN
        READ(10,*) dimx, dimy, dimz
      ELSE IF (key3 == 'int' .OR. key3 == 'INT') THEN
        READ(10,*) inter
      ELSE IF (key3 == 'cro' .OR. key3 == 'CRO') THEN
        READ(10,*) cross
      ELSE IF (key3 == 'apr' .OR. key3 == 'APR') THEN
        READ(10,*) aproxima
      ELSE IF (key3 == 'end' .OR. key3 == 'END') THEN
        EXIT
      ENDIF
    ENDDO
    CLOSE(10)

  END SUBROUTINE read_ini

END PROGRAM prepara
