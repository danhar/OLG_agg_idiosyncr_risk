module ranking_full_mod
! ranks an array
! see http://www.fortran-2000.com/rank/ , Subroutine MRGRNK (XVALT, IMULT)
    implicit none
    integer, parameter :: kdp = selected_real_kind(15)
    public  :: rank
    private :: kdp
    private :: R_mrgrnk, I_mrgrnk, D_mrgrnk

    interface rank
      module procedure D_mrgrnk, R_mrgrnk, I_mrgrnk
    end interface rank

contains

pure subroutine D_mrgrnk (XDONT, IRNGT)
! __________________________________________________________
!   MRGRNK = Merge-sort ranking of an array
!   For performance reasons, the first 2 passes are taken
!   out of the standard loop, and use dedicated coding.
! __________________________________________________________
! __________________________________________________________
      real (kind=kdp), dimension (:), intent (in) :: XDONT
      integer, dimension (:), intent (out) :: IRNGT
! __________________________________________________________
      real (kind=kdp) :: XVALA, XVALB
!
      integer, dimension (SIZE(IRNGT)) :: JWRKT
      integer :: LMTNA, LMTNC, IRNG1, IRNG2
      integer :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
!
      NVAL = Min (SIZE(XDONT), SIZE(IRNGT))
      select case (NVAL)
      case (:0)
         return
      case (1)
         IRNGT (1) = 1
         return
      case default
         continue
      end select
!
!  Fill-in the index array, creating ordered couples
!
      do IIND = 2, NVAL, 2
         if (XDONT(IIND-1) <= XDONT(IIND)) then
            IRNGT (IIND-1) = IIND - 1
            IRNGT (IIND) = IIND
         else
            IRNGT (IIND-1) = IIND
            IRNGT (IIND) = IIND - 1
         end if
      end do
      if (Modulo(NVAL, 2) /= 0) then
         IRNGT (NVAL) = NVAL
      end if
!
!  We will now have ordered subsets A - B - A - B - ...
!  and merge A and B couples into     C   -   C   - ...
!
      LMTNA = 2
      LMTNC = 4
!
!  First iteration. The length of the ordered subsets goes from 2 to 4
!
      do
         if (NVAL <= 2) exit
!
!   Loop on merges of A and B into C
!
         do IWRKD = 0, NVAL - 1, 4
            if ((IWRKD+4) > NVAL) then
               if ((IWRKD+2) >= NVAL) exit
!
!   1 2 3
!
               if (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) exit
!
!   1 3 2
!
               if (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) then
                  IRNG2 = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNG2
!
!   3 1 2
!
               else
                  IRNG1 = IRNGT (IWRKD+1)
                  IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNG1
               end if
               exit
            end if
!
!   1 2 3 4
!
            if (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) cycle
!
!   1 3 x x
!
            if (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) then
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
               if (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) then
!   1 3 2 4
                  IRNGT (IWRKD+3) = IRNG2
               else
!   1 3 4 2
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+4) = IRNG2
               end if
!
!   3 x x x
!
            else
               IRNG1 = IRNGT (IWRKD+1)
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
               if (XDONT(IRNG1) <= XDONT(IRNGT(IWRKD+4))) then
                  IRNGT (IWRKD+2) = IRNG1
                  if (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) then
!   3 1 2 4
                     IRNGT (IWRKD+3) = IRNG2
                  else
!   3 1 4 2
                     IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                     IRNGT (IWRKD+4) = IRNG2
                  end if
               else
!   3 4 1 2
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+3) = IRNG1
                  IRNGT (IWRKD+4) = IRNG2
               end if
            end if
         end do
!
!  The Cs become As and Bs
!
         LMTNA = 4
         exit
      end do
!
!  Iteration loop. Each time, the length of the ordered subsets
!  is doubled.
!
      do
         if (LMTNA >= NVAL) exit
         IWRKF = 0
         LMTNC = 2 * LMTNC
!
!   Loop on merges of A and B into C
!
         do
            IWRK = IWRKF
            IWRKD = IWRKF + 1
            JINDA = IWRKF + LMTNA
            IWRKF = IWRKF + LMTNC
            if (IWRKF >= NVAL) then
               if (JINDA >= NVAL) exit
               IWRKF = NVAL
            end if
            IINDA = 1
            IINDB = JINDA + 1
!
!   Shortcut for the case when the max of A is smaller
!   than the min of B. This line may be activated when the
!   initial set is already close to sorted.
!
!          IF (XDONT(IRNGT(JINDA)) <= XDONT(IRNGT(IINDB))) CYCLE
!
!  One steps in the C subset, that we build in the final rank array
!
!  Make a copy of the rank array for the merge iteration
!
            JWRKT (1:LMTNA) = IRNGT (IWRKD:JINDA)
!
            XVALA = XDONT (JWRKT(IINDA))
            XVALB = XDONT (IRNGT(IINDB))
!
            do
               IWRK = IWRK + 1
!
!  We still have unprocessed values in both A and B
!
               if (XVALA > XVALB) then
                  IRNGT (IWRK) = IRNGT (IINDB)
                  IINDB = IINDB + 1
                  if (IINDB > IWRKF) then
!  Only A still with unprocessed values
                     IRNGT (IWRK+1:IWRKF) = JWRKT (IINDA:LMTNA)
                     exit
                  end if
                  XVALB = XDONT (IRNGT(IINDB))
               else
                  IRNGT (IWRK) = JWRKT (IINDA)
                  IINDA = IINDA + 1
                  if (IINDA > LMTNA) exit! Only B still with unprocessed values
                  XVALA = XDONT (JWRKT(IINDA))
               end if
!
            end do
         end do
!
!  The Cs become As and Bs
!
         LMTNA = 2 * LMTNA
      end do
!
      return
!
end subroutine D_mrgrnk

pure subroutine R_mrgrnk (XDONT, IRNGT)
! __________________________________________________________
!   MRGRNK = Merge-sort ranking of an array
!   For performance reasons, the first 2 passes are taken
!   out of the standard loop, and use dedicated coding.
! __________________________________________________________
! _________________________________________________________
      real, dimension (:), intent (in) :: XDONT
      integer, dimension (:), intent (out) :: IRNGT
! __________________________________________________________
      real :: XVALA, XVALB
!
      integer, dimension (SIZE(IRNGT)) :: JWRKT
      integer :: LMTNA, LMTNC, IRNG1, IRNG2
      integer :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
!
      NVAL = Min (SIZE(XDONT), SIZE(IRNGT))
      select case (NVAL)
      case (:0)
         return
      case (1)
         IRNGT (1) = 1
         return
      case default
         continue
      end select
!
!  Fill-in the index array, creating ordered couples
!
      do IIND = 2, NVAL, 2
         if (XDONT(IIND-1) <= XDONT(IIND)) then
            IRNGT (IIND-1) = IIND - 1
            IRNGT (IIND) = IIND
         else
            IRNGT (IIND-1) = IIND
            IRNGT (IIND) = IIND - 1
         end if
      end do
      if (Modulo(NVAL, 2) /= 0) then
         IRNGT (NVAL) = NVAL
      end if
!
!  We will now have ordered subsets A - B - A - B - ...
!  and merge A and B couples into     C   -   C   - ...
!
      LMTNA = 2
      LMTNC = 4
!
!  First iteration. The length of the ordered subsets goes from 2 to 4
!
      do
         if (NVAL <= 2) exit
!
!   Loop on merges of A and B into C
!
         do IWRKD = 0, NVAL - 1, 4
            if ((IWRKD+4) > NVAL) then
               if ((IWRKD+2) >= NVAL) exit
!
!   1 2 3
!
               if (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) exit
!
!   1 3 2
!
               if (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) then
                  IRNG2 = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNG2
!
!   3 1 2
!
               else
                  IRNG1 = IRNGT (IWRKD+1)
                  IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNG1
               end if
               exit
            end if
!
!   1 2 3 4
!
            if (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) cycle
!
!   1 3 x x
!
            if (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) then
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
               if (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) then
!   1 3 2 4
                  IRNGT (IWRKD+3) = IRNG2
               else
!   1 3 4 2
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+4) = IRNG2
               end if
!
!   3 x x x
!
            else
               IRNG1 = IRNGT (IWRKD+1)
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
               if (XDONT(IRNG1) <= XDONT(IRNGT(IWRKD+4))) then
                  IRNGT (IWRKD+2) = IRNG1
                  if (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) then
!   3 1 2 4
                     IRNGT (IWRKD+3) = IRNG2
                  else
!   3 1 4 2
                     IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                     IRNGT (IWRKD+4) = IRNG2
                  end if
               else
!   3 4 1 2
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+3) = IRNG1
                  IRNGT (IWRKD+4) = IRNG2
               end if
            end if
         end do
!
!  The Cs become As and Bs
!
         LMTNA = 4
         exit
      end do
!
!  Iteration loop. Each time, the length of the ordered subsets
!  is doubled.
!
      do
         if (LMTNA >= NVAL) exit
         IWRKF = 0
         LMTNC = 2 * LMTNC
!
!   Loop on merges of A and B into C
!
         do
            IWRK = IWRKF
            IWRKD = IWRKF + 1
            JINDA = IWRKF + LMTNA
            IWRKF = IWRKF + LMTNC
            if (IWRKF >= NVAL) then
               if (JINDA >= NVAL) exit
               IWRKF = NVAL
            end if
            IINDA = 1
            IINDB = JINDA + 1
!
!   Shortcut for the case when the max of A is smaller
!   than the min of B. This line may be activated when the
!   initial set is already close to sorted.
!
!          IF (XDONT(IRNGT(JINDA)) <= XDONT(IRNGT(IINDB))) CYCLE
!
!  One steps in the C subset, that we build in the final rank array
!
!  Make a copy of the rank array for the merge iteration
!
            JWRKT (1:LMTNA) = IRNGT (IWRKD:JINDA)
!
            XVALA = XDONT (JWRKT(IINDA))
            XVALB = XDONT (IRNGT(IINDB))
!
            do
               IWRK = IWRK + 1
!
!  We still have unprocessed values in both A and B
!
               if (XVALA > XVALB) then
                  IRNGT (IWRK) = IRNGT (IINDB)
                  IINDB = IINDB + 1
                  if (IINDB > IWRKF) then
!  Only A still with unprocessed values
                     IRNGT (IWRK+1:IWRKF) = JWRKT (IINDA:LMTNA)
                     exit
                  end if
                  XVALB = XDONT (IRNGT(IINDB))
               else
                  IRNGT (IWRK) = JWRKT (IINDA)
                  IINDA = IINDA + 1
                  if (IINDA > LMTNA) exit! Only B still with unprocessed values
                  XVALA = XDONT (JWRKT(IINDA))
               end if
!
            end do
         end do
!
!  The Cs become As and Bs
!
         LMTNA = 2 * LMTNA
      end do
!
      return
!
end subroutine R_mrgrnk

pure subroutine I_mrgrnk (XDONT, IRNGT)
! __________________________________________________________
!   MRGRNK = Merge-sort ranking of an array
!   For performance reasons, the first 2 passes are taken
!   out of the standard loop, and use dedicated coding.
! __________________________________________________________
! __________________________________________________________
      integer, dimension (:), intent (in)  :: XDONT
      integer, dimension (:), intent (out) :: IRNGT
! __________________________________________________________
      integer :: XVALA, XVALB
!
      integer, dimension (SIZE(IRNGT)) :: JWRKT
      integer :: LMTNA, LMTNC, IRNG1, IRNG2
      integer :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
!
      NVAL = Min (SIZE(XDONT), SIZE(IRNGT))
      select case (NVAL)
      case (:0)
         return
      case (1)
         IRNGT (1) = 1
         return
      case default
         continue
      end select
!
!  Fill-in the index array, creating ordered couples
!
      do IIND = 2, NVAL, 2
         if (XDONT(IIND-1) <= XDONT(IIND)) then
            IRNGT (IIND-1) = IIND - 1
            IRNGT (IIND) = IIND
         else
            IRNGT (IIND-1) = IIND
            IRNGT (IIND) = IIND - 1
         end if
      end do
      if (Modulo(NVAL, 2) /= 0) then
         IRNGT (NVAL) = NVAL
      end if
!
!  We will now have ordered subsets A - B - A - B - ...
!  and merge A and B couples into     C   -   C   - ...
!
      LMTNA = 2
      LMTNC = 4
!
!  First iteration. The length of the ordered subsets goes from 2 to 4
!
      do
         if (NVAL <= 2) exit
!
!   Loop on merges of A and B into C
!
         do IWRKD = 0, NVAL - 1, 4
            if ((IWRKD+4) > NVAL) then
               if ((IWRKD+2) >= NVAL) exit
!
!   1 2 3
!
               if (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) exit
!
!   1 3 2
!
               if (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) then
                  IRNG2 = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNG2
!
!   3 1 2
!
               else
                  IRNG1 = IRNGT (IWRKD+1)
                  IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNG1
               end if
               exit
            end if
!
!   1 2 3 4
!
            if (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) cycle
!
!   1 3 x x
!
            if (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) then
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
               if (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) then
!   1 3 2 4
                  IRNGT (IWRKD+3) = IRNG2
               else
!   1 3 4 2
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+4) = IRNG2
               end if
!
!   3 x x x
!
            else
               IRNG1 = IRNGT (IWRKD+1)
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
               if (XDONT(IRNG1) <= XDONT(IRNGT(IWRKD+4))) then
                  IRNGT (IWRKD+2) = IRNG1
                  if (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) then
!   3 1 2 4
                     IRNGT (IWRKD+3) = IRNG2
                  else
!   3 1 4 2
                     IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                     IRNGT (IWRKD+4) = IRNG2
                  end if
               else
!   3 4 1 2
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+3) = IRNG1
                  IRNGT (IWRKD+4) = IRNG2
               end if
            end if
         end do
!
!  The Cs become As and Bs
!
         LMTNA = 4
         exit
      end do
!
!  Iteration loop. Each time, the length of the ordered subsets
!  is doubled.
!
      do
         if (LMTNA >= NVAL) exit
         IWRKF = 0
         LMTNC = 2 * LMTNC
!
!   Loop on merges of A and B into C
!
         do
            IWRK = IWRKF
            IWRKD = IWRKF + 1
            JINDA = IWRKF + LMTNA
            IWRKF = IWRKF + LMTNC
            if (IWRKF >= NVAL) then
               if (JINDA >= NVAL) exit
               IWRKF = NVAL
            end if
            IINDA = 1
            IINDB = JINDA + 1
!
!   Shortcut for the case when the max of A is smaller
!   than the min of B. This line may be activated when the
!   initial set is already close to sorted.
!
!          IF (XDONT(IRNGT(JINDA)) <= XDONT(IRNGT(IINDB))) CYCLE
!
!  One steps in the C subset, that we build in the final rank array
!
!  Make a copy of the rank array for the merge iteration
!
            JWRKT (1:LMTNA) = IRNGT (IWRKD:JINDA)
!
            XVALA = XDONT (JWRKT(IINDA))
            XVALB = XDONT (IRNGT(IINDB))
!
            do
               IWRK = IWRK + 1
!
!  We still have unprocessed values in both A and B
!
               if (XVALA > XVALB) then
                  IRNGT (IWRK) = IRNGT (IINDB)
                  IINDB = IINDB + 1
                  if (IINDB > IWRKF) then
!  Only A still with unprocessed values
                     IRNGT (IWRK+1:IWRKF) = JWRKT (IINDA:LMTNA)
                     exit
                  end if
                  XVALB = XDONT (IRNGT(IINDB))
               else
                  IRNGT (IWRK) = JWRKT (IINDA)
                  IINDA = IINDA + 1
                  if (IINDA > LMTNA) exit! Only B still with unprocessed values
                  XVALA = XDONT (JWRKT(IINDA))
               end if
!
            end do
         end do
!
!  The Cs become As and Bs
!
         LMTNA = 2 * LMTNA
      end do
!
      return
!
end subroutine I_mrgrnk

end module ranking_full_mod
