module sorting_full_mod
! see http://www.fortran-2000.com/rank/ , Subroutine REFSOR (XVALT)
! sorts an array
    use kinds ,only : kdp => dp

    implicit none
    private

    public  :: sort
    private :: R_refsor, I_refsor, D_refsor
    private :: R_inssor, I_inssor, D_inssor
    private :: R_subsor, I_subsor, D_subsor

    interface sort
      module procedure d_refsor, r_refsor, i_refsor
    end interface sort

contains

pure subroutine D_refsor (XDONT)
!  Sorts XDONT into ascending order - Quicksort
! __________________________________________________________
!  Quicksort chooses a "pivot" in the set, and explores the
!  array from both ends, looking for a value > pivot with the
!  increasing index, for a value <= pivot with the decreasing
!  index, and swapping them when it has found one of each.
!  The array is then subdivided in 2 ([3]) subsets:
!  { values <= pivot} {pivot} {values > pivot}
!  One then call recursively the program to sort each subset.
!  When the size of the subarray is small enough, one uses an
!  insertion sort that is faster for very small sets.
!  Michel Olagnon - Apr. 2000
! __________________________________________________________
! __________________________________________________________
      real (kind=kdp), dimension (:), intent (inout) :: XDONT
! __________________________________________________________
!
!
      call D_subsor (XDONT, 1, Size (XDONT))
      call D_inssor (XDONT)
      return
end subroutine D_refsor

pure recursive subroutine D_subsor (XDONT, IDEB1, IFIN1)
!  Sorts XDONT from IDEB1 to IFIN1
! __________________________________________________________
      real(kind=kdp), dimension (:), intent (inout) :: XDONT
      integer, intent (in) :: IDEB1, IFIN1
! __________________________________________________________
      integer, parameter :: NINS = 16 ! Max for insertion sort
      integer :: ICRS, IDEB, IDCR, IFIN, IMIL
      real(kind=kdp) :: XPIV, XWRK
!
      IDEB = IDEB1
      IFIN = IFIN1
!
!  If we don't have enough values to make it worth while, we leave
!  them unsorted, and the final insertion sort will take care of them
!
      if ((IFIN - IDEB) > NINS) then
         IMIL = (IDEB+IFIN) / 2
!
!  One chooses a pivot, median of 1st, last, and middle values
!
         if (XDONT(IMIL) < XDONT(IDEB)) then
            XWRK = XDONT (IDEB)
            XDONT (IDEB) = XDONT (IMIL)
            XDONT (IMIL) = XWRK
         end if
         if (XDONT(IMIL) > XDONT(IFIN)) then
            XWRK = XDONT (IFIN)
            XDONT (IFIN) = XDONT (IMIL)
            XDONT (IMIL) = XWRK
            if (XDONT(IMIL) < XDONT(IDEB)) then
               XWRK = XDONT (IDEB)
               XDONT (IDEB) = XDONT (IMIL)
               XDONT (IMIL) = XWRK
            end if
         end if
         XPIV = XDONT (IMIL)
!
!  One exchanges values to put those > pivot in the end and
!  those <= pivot at the beginning
!
         ICRS = IDEB
         IDCR = IFIN
         ECH2: do
            do
               ICRS = ICRS + 1
               if (ICRS >= IDCR) then
!
!  the first  >  pivot is IDCR
!  the last   <= pivot is ICRS-1
!  Note: If one arrives here on the first iteration, then
!        the pivot is the maximum of the set, the last value is equal
!        to it, and one can reduce by one the size of the set to process,
!        as if XDONT (IFIN) > XPIV
!
                  exit ECH2
!
               end if
               if (XDONT(ICRS) > XPIV) exit
            end do
            do
               if (XDONT(IDCR) <= XPIV) exit
               IDCR = IDCR - 1
               if (ICRS >= IDCR) then
!
!  The last value < pivot is always ICRS-1
!
                  exit ECH2
               end if
            end do
!
            XWRK = XDONT (IDCR)
            XDONT (IDCR) = XDONT (ICRS)
            XDONT (ICRS) = XWRK
         end do ECH2
!
!  One now sorts each of the two sub-intervals
!
         call D_subsor (XDONT, IDEB1, ICRS-1)
         call D_subsor (XDONT, IDCR, IFIN1)
      end if
      return
   end subroutine D_subsor

   pure subroutine D_inssor (XDONT)
!  Sorts XDONT into increasing order (Insertion sort)
! __________________________________________________________
      real(kind=kdp), dimension (:), intent (inout) :: XDONT
! __________________________________________________________
      integer :: ICRS, IDCR
      real(kind=kdp) :: XWRK
!
      do ICRS = 2, Size (XDONT)
         XWRK = XDONT (ICRS)
         if (XWRK >= XDONT(ICRS-1)) cycle
         XDONT (ICRS) = XDONT (ICRS-1)
         do IDCR = ICRS - 2, 1, - 1
            if (XWRK >= XDONT(IDCR)) exit
            XDONT (IDCR+1) = XDONT (IDCR)
         end do
         XDONT (IDCR+1) = XWRK
      end do
!
      return
!
end subroutine D_inssor
!
pure subroutine R_refsor (XDONT)
!  Sorts XDONT into ascending order - Quicksort
! __________________________________________________________
!  Quicksort chooses a "pivot" in the set, and explores the
!  array from both ends, looking for a value > pivot with the
!  increasing index, for a value <= pivot with the decreasing
!  index, and swapping them when it has found one of each.
!  The array is then subdivided in 2 ([3]) subsets:
!  { values <= pivot} {pivot} {values > pivot}
!  One then call recursively the program to sort each subset.
!  When the size of the subarray is small enough, one uses an
!  insertion sort that is faster for very small sets.
!  Michel Olagnon - Apr. 2000
! __________________________________________________________
! _________________________________________________________
      real, dimension (:), intent (inout) :: XDONT
! __________________________________________________________
!
!
      call R_subsor (XDONT, 1, Size (XDONT))
      call R_inssor (XDONT)
      return
end subroutine R_refsor
pure recursive subroutine R_subsor (XDONT, IDEB1, IFIN1)
!  Sorts XDONT from IDEB1 to IFIN1
! __________________________________________________________
      real, dimension (:), intent (inout) :: XDONT
      integer, intent (in) :: IDEB1, IFIN1
! __________________________________________________________
      integer, parameter :: NINS = 16 ! Max for insertion sort
      integer :: ICRS, IDEB, IDCR, IFIN, IMIL
      real :: XPIV, XWRK
!
      IDEB = IDEB1
      IFIN = IFIN1
!
!  If we don't have enough values to make it worth while, we leave
!  them unsorted, and the final insertion sort will take care of them
!
      if ((IFIN - IDEB) > NINS) then
         IMIL = (IDEB+IFIN) / 2
!
!  One chooses a pivot, median of 1st, last, and middle values
!
         if (XDONT(IMIL) < XDONT(IDEB)) then
            XWRK = XDONT (IDEB)
            XDONT (IDEB) = XDONT (IMIL)
            XDONT (IMIL) = XWRK
         end if
         if (XDONT(IMIL) > XDONT(IFIN)) then
            XWRK = XDONT (IFIN)
            XDONT (IFIN) = XDONT (IMIL)
            XDONT (IMIL) = XWRK
            if (XDONT(IMIL) < XDONT(IDEB)) then
               XWRK = XDONT (IDEB)
               XDONT (IDEB) = XDONT (IMIL)
               XDONT (IMIL) = XWRK
            end if
         end if
         XPIV = XDONT (IMIL)
!
!  One exchanges values to put those > pivot in the end and
!  those <= pivot at the beginning
!
         ICRS = IDEB
         IDCR = IFIN
         ECH2: do
            do
               ICRS = ICRS + 1
               if (ICRS >= IDCR) then
!
!  the first  >  pivot is IDCR
!  the last   <= pivot is ICRS-1
!  Note: If one arrives here on the first iteration, then
!        the pivot is the maximum of the set, the last value is equal
!        to it, and one can reduce by one the size of the set to process,
!        as if XDONT (IFIN) > XPIV
!
                  exit ECH2
!
               end if
               if (XDONT(ICRS) > XPIV) exit
            end do
            do
               if (XDONT(IDCR) <= XPIV) exit
               IDCR = IDCR - 1
               if (ICRS >= IDCR) then
!
!  The last value < pivot is always ICRS-1
!
                  exit ECH2
               end if
            end do
!
            XWRK = XDONT (IDCR)
            XDONT (IDCR) = XDONT (ICRS)
            XDONT (ICRS) = XWRK
         end do ECH2
!
!  One now sorts each of the two sub-intervals
!
         call R_subsor (XDONT, IDEB1, ICRS-1)
         call R_subsor (XDONT, IDCR, IFIN1)
      end if
      return
   end subroutine R_subsor
   pure subroutine R_inssor (XDONT)
!  Sorts XDONT into increasing order (Insertion sort)
! __________________________________________________________
      real, dimension (:), intent (inout) :: XDONT
! __________________________________________________________
      integer :: ICRS, IDCR
      real :: XWRK
!
      do ICRS = 2, Size (XDONT)
         XWRK = XDONT (ICRS)
         if (XWRK >= XDONT(ICRS-1)) cycle
         XDONT (ICRS) = XDONT (ICRS-1)
         do IDCR = ICRS - 2, 1, - 1
            if (XWRK >= XDONT(IDCR)) exit
            XDONT (IDCR+1) = XDONT (IDCR)
         end do
         XDONT (IDCR+1) = XWRK
      end do
!
      return
!
end subroutine R_inssor
!
pure subroutine I_refsor (XDONT)
!  Sorts XDONT into ascending order - Quicksort
! __________________________________________________________
!  Quicksort chooses a "pivot" in the set, and explores the
!  array from both ends, looking for a value > pivot with the
!  increasing index, for a value <= pivot with the decreasing
!  index, and swapping them when it has found one of each.
!  The array is then subdivided in 2 ([3]) subsets:
!  { values <= pivot} {pivot} {values > pivot}
!  One then call recursively the program to sort each subset.
!  When the size of the subarray is small enough, one uses an
!  insertion sort that is faster for very small sets.
!  Michel Olagnon - Apr. 2000
! __________________________________________________________
! __________________________________________________________
      integer, dimension (:), intent (inout)  :: XDONT
! __________________________________________________________
!
!
      call I_subsor (XDONT, 1, Size (XDONT))
      call I_inssor (XDONT)
      return
end subroutine I_refsor

pure recursive subroutine I_subsor (XDONT, IDEB1, IFIN1)
!  Sorts XDONT from IDEB1 to IFIN1
! __________________________________________________________
      integer, dimension (:), intent (inout) :: XDONT
      integer, intent (in) :: IDEB1, IFIN1
! __________________________________________________________
      integer, parameter :: NINS = 16 ! Max for insertion sort
      integer :: ICRS, IDEB, IDCR, IFIN, IMIL
      integer :: XPIV, XWRK
!
      IDEB = IDEB1
      IFIN = IFIN1
!
!  If we don't have enough values to make it worth while, we leave
!  them unsorted, and the final insertion sort will take care of them
!
      if ((IFIN - IDEB) > NINS) then
         IMIL = (IDEB+IFIN) / 2
!
!  One chooses a pivot, median of 1st, last, and middle values
!
         if (XDONT(IMIL) < XDONT(IDEB)) then
            XWRK = XDONT (IDEB)
            XDONT (IDEB) = XDONT (IMIL)
            XDONT (IMIL) = XWRK
         end if
         if (XDONT(IMIL) > XDONT(IFIN)) then
            XWRK = XDONT (IFIN)
            XDONT (IFIN) = XDONT (IMIL)
            XDONT (IMIL) = XWRK
            if (XDONT(IMIL) < XDONT(IDEB)) then
               XWRK = XDONT (IDEB)
               XDONT (IDEB) = XDONT (IMIL)
               XDONT (IMIL) = XWRK
            end if
         end if
         XPIV = XDONT (IMIL)
!
!  One exchanges values to put those > pivot in the end and
!  those <= pivot at the beginning
!
         ICRS = IDEB
         IDCR = IFIN
         ECH2: do
            do
               ICRS = ICRS + 1
               if (ICRS >= IDCR) then
!
!  the first  >  pivot is IDCR
!  the last   <= pivot is ICRS-1
!  Note: If one arrives here on the first iteration, then
!        the pivot is the maximum of the set, the last value is equal
!        to it, and one can reduce by one the size of the set to process,
!        as if XDONT (IFIN) > XPIV
!
                  exit ECH2
!
               end if
               if (XDONT(ICRS) > XPIV) exit
            end do
            do
               if (XDONT(IDCR) <= XPIV) exit
               IDCR = IDCR - 1
               if (ICRS >= IDCR) then
!
!  The last value < pivot is always ICRS-1
!
                  exit ECH2
               end if
            end do
!
            XWRK = XDONT (IDCR)
            XDONT (IDCR) = XDONT (ICRS)
            XDONT (ICRS) = XWRK
         end do ECH2
!
!  One now sorts each of the two sub-intervals
!
         call I_subsor (XDONT, IDEB1, ICRS-1)
         call I_subsor (XDONT, IDCR, IFIN1)
      end if
      return
   end subroutine I_subsor

   pure subroutine I_inssor (XDONT)
!  Sorts XDONT into increasing order (Insertion sort)
! __________________________________________________________
      integer, dimension (:), intent (inout) :: XDONT
! __________________________________________________________
      integer :: ICRS, IDCR
      integer :: XWRK
!
      do ICRS = 2, Size (XDONT)
         XWRK = XDONT (ICRS)
         if (XWRK >= XDONT(ICRS-1)) cycle
         XDONT (ICRS) = XDONT (ICRS-1)
         do IDCR = ICRS - 2, 1, - 1
            if (XWRK >= XDONT(IDCR)) exit
            XDONT (IDCR+1) = XDONT (IDCR)
         end do
         XDONT (IDCR+1) = XWRK
      end do
!
      return
!
end subroutine I_inssor
!
end module sorting_full_mod
