module sorting_partial_mod
! generalizes concept of Median
! see http://www.fortran-2000.com/rank/
    use kinds ,only : kdp => dp
    implicit none

    public :: valnth, indnth
    private :: R_valnth, I_valnth, D_valnth

    interface valnth
      module procedure d_valnth, r_valnth, i_valnth
    end interface valnth

    interface indnth
      module procedure d_indnth, r_indnth, i_indnth
    end interface indnth

contains

pure function D_valnth (XDONT, NORD) result (valnth)
!  Return NORDth value of XDONT, i.e fractile of order NORD/SIZE(XDONT).
! __________________________________________________________
!  This routine uses a pivoting strategy such as the one of
!  finding the median based on the quicksort algorithm, but
!  we skew the pivot choice to try to bring it to NORD as
!  fast as possible. It uses 2 temporary arrays, where it
!  stores the indices of the values smaller than the pivot
!  (ILOWT), and the indices of values larger than the pivot
!  that we might still need later on (IHIGT). It iterates
!  until it can bring the number of values in ILOWT to
!  exactly NORD, and then finds the maximum of this set.
!  Michel Olagnon - Aug. 2000
! __________________________________________________________
! __________________________________________________________
      real (kind=kdp), dimension (:), intent (in) :: XDONT
      real (kind=kdp) :: valnth
      integer, intent (in) :: NORD
! __________________________________________________________
      real (kind=kdp), dimension (SIZE(XDONT)) :: XLOWT, XHIGT
      real (kind=kdp) :: XPIV, XPIV0, XWRK, XWRK1, XWRK2, XWRK3, XMIN, XMAX
!
      integer :: NDON, JHIG, JLOW, IHIG
      integer :: IMIL, IFIN, ICRS, IDCR, ILOW
      integer :: JLM2, JLM1, JHM2, JHM1, INTH
!
      NDON = SIZE (XDONT)
      INTH = MAX (MIN (NORD, NDON), 1)
!
!    First loop is used to fill-in XLOWT, XHIGT at the same time
!
      if (NDON < 2) then
         if (INTH == 1) VALNTH = XDONT (1)
         return
      end if
!
!  One chooses a pivot, best estimate possible to put fractile near
!  mid-point of the set of low values.
!
      if (XDONT(2) < XDONT(1)) then
         XLOWT (1) = XDONT(2)
         XHIGT (1) = XDONT(1)
      else
         XLOWT (1) = XDONT(1)
         XHIGT (1) = XDONT(2)
      end if
!
      if (NDON < 3) then
         if (INTH == 1) VALNTH = XLOWT (1)
         if (INTH == 2) VALNTH = XHIGT (1)
         return
      end if
!
      if (XDONT(3) < XHIGT(1)) then
         XHIGT (2) = XHIGT (1)
         if (XDONT(3) < XLOWT(1)) then
            XHIGT (1) = XLOWT (1)
            XLOWT (1) = XDONT(3)
         else
            XHIGT (1) = XDONT(3)
         end if
      else
         XHIGT (2) = XDONT(3)
      end if
!
      if (NDON < 4) then
         if (INTH == 1) then
             VALNTH = XLOWT (1)
         else
             VALNTH = XHIGT (INTH - 1)
         end if
         return
      end if
!
      if (XDONT(NDON) < XHIGT(1)) then
         XHIGT (3) = XHIGT (2)
         XHIGT (2) = XHIGT (1)
         if (XDONT(NDON) < XLOWT(1)) then
            XHIGT (1) = XLOWT (1)
            XLOWT (1) = XDONT(NDON)
         else
            XHIGT (1) = XDONT(NDON)
         end if
      else
         XHIGT (3) = XDONT(NDON)
      end if
!
      if (NDON < 5) then
         if (INTH == 1) then
             VALNTH = XLOWT (1)
         else
             VALNTH = XHIGT (INTH - 1)
         end if
         return
      end if
!

      JLOW = 1
      JHIG = 3
      XPIV = XLOWT(1) + REAL(2*INTH)/REAL(NDON+INTH) * (XHIGT(3)-XLOWT(1))
      if (XPIV >= XHIGT(1)) then
         XPIV = XLOWT(1) + REAL(2*INTH)/REAL(NDON+INTH) * &
                           (XHIGT(2)-XLOWT(1))
         if (XPIV >= XHIGT(1)) &
             XPIV = XLOWT(1) + REAL (2*INTH) / REAL (NDON+INTH) * &
                               (XHIGT(1)-XLOWT(1))
      end if
      XPIV0 = XPIV
!
!  One puts values > pivot in the end and those <= pivot
!  at the beginning. This is split in 2 cases, so that
!  we can skip the loop test a number of times.
!  As we are also filling in the work arrays at the same time
!  we stop filling in the XHIGT array as soon as we have more
!  than enough values in XLOWT.
!
!
      if (XDONT(NDON) > XPIV) then
         ICRS = 3
         do
            ICRS = ICRS + 1
            if (XDONT(ICRS) > XPIV) then
               if (ICRS >= NDON) exit
               JHIG = JHIG + 1
               XHIGT (JHIG) = XDONT(ICRS)
            else
               JLOW = JLOW + 1
               XLOWT (JLOW) = XDONT(ICRS)
               if (JLOW >= INTH) exit
            end if
         end do
!
!  One restricts further processing because it is no use
!  to store more high values
!
         if (ICRS < NDON-1) then
            do
               ICRS = ICRS + 1
               if (XDONT(ICRS) <= XPIV) then
                  JLOW = JLOW + 1
                  XLOWT (JLOW) = XDONT(ICRS)
               else if (ICRS >= NDON) then
                  exit
               end if
            end do
         end if
!
!
      else
!
!  Same as above, but this is not as easy to optimize, so the
!  DO-loop is kept
!
         do ICRS = 4, NDON - 1
            if (XDONT(ICRS) > XPIV) then
               JHIG = JHIG + 1
               XHIGT (JHIG) = XDONT(ICRS)
            else
               JLOW = JLOW + 1
               XLOWT (JLOW) = XDONT(ICRS)
               if (JLOW >= INTH) exit
            end if
         end do
!
         if (ICRS < NDON-1) then
            do
               ICRS = ICRS + 1
               if (XDONT(ICRS) <= XPIV) then
                  if (ICRS >= NDON) exit
                  JLOW = JLOW + 1
                  XLOWT (JLOW) = XDONT(ICRS)
               end if
            end do
         end if
      end if
!
      JLM2 = 0
      JLM1 = 0
      JHM2 = 0
      JHM1 = 0
      do
         if (JLM2 == JLOW .and. JHM2 == JHIG) then
!
!   We are oscillating. Perturbate by bringing JLOW closer by one
!   to INTH
!
             if (INTH > JLOW) then
                XMIN = XHIGT(1)
                IHIG = 1
                do ICRS = 2, JHIG
                   if (XHIGT(ICRS) < XMIN) then
                      XMIN = XHIGT(ICRS)
                      IHIG = ICRS
                   end if
                end do
!
                JLOW = JLOW + 1
                XLOWT (JLOW) = XHIGT (IHIG)
                XHIGT (IHIG) = XHIGT (JHIG)
                JHIG = JHIG - 1
             else

                XMAX = XLOWT (JLOW)
                JLOW = JLOW - 1
                do ICRS = 1, JLOW
                   if (XLOWT(ICRS) > XMAX) then
                      XWRK = XMAX
                      XMAX = XLOWT(ICRS)
                      XLOWT (ICRS) = XWRK
                   end if
                end do
             end if
         end if
         JLM2 = JLM1
         JLM1 = JLOW
         JHM2 = JHM1
         JHM1 = JHIG
!
!   We try to bring the number of values in the low values set
!   closer to INTH.
!
         select case (INTH-JLOW)
         case (2:)
!
!   Not enough values in low part, at least 2 are missing
!
            INTH = INTH - JLOW
            JLOW = 0
            select case (JHIG)
!!!!!           CASE DEFAULT
!!!!!              write (unit=*,fmt=*) "Assertion failed"
!!!!!              STOP
!
!   We make a special case when we have so few values in
!   the high values set that it is bad performance to choose a pivot
!   and apply the general algorithm.
!
            case (2)
               if (XHIGT(1) <= XHIGT(2)) then
                  JLOW = JLOW + 1
                  XLOWT (JLOW) = XHIGT (1)
                  JLOW = JLOW + 1
                  XLOWT (JLOW) = XHIGT (2)
               else
                  JLOW = JLOW + 1
                  XLOWT (JLOW) = XHIGT (2)
                  JLOW = JLOW + 1
                  XLOWT (JLOW) = XHIGT (1)
               end if
               exit
!
            case (3)
!
!
               XWRK1 = XHIGT (1)
               XWRK2 = XHIGT (2)
               XWRK3 = XHIGT (3)
               if (XWRK2 < XWRK1) then
                  XHIGT (1) = XWRK2
                  XHIGT (2) = XWRK1
                  XWRK2 = XWRK1
               end if
               if (XWRK2 > XWRK3) then
                  XHIGT (3) = XWRK2
                  XHIGT (2) = XWRK3
                  XWRK2 = XWRK3
                  if (XWRK2 < XHIGT(1)) then
                     XHIGT (2) = XHIGT (1)
                     XHIGT (1) = XWRK2
                  end if
               end if
               JHIG = 0
               do ICRS = JLOW + 1, INTH
                  JHIG = JHIG + 1
                  XLOWT (ICRS) = XHIGT (JHIG)
               end do
               JLOW = INTH
               exit
!
            case (4:)
!
!
               XPIV0 = XPIV
               IFIN = JHIG
!
!  One chooses a pivot from the 2 first values and the last one.
!  This should ensure sufficient renewal between iterations to
!  avoid worst case behavior effects.
!
               XWRK1 = XHIGT (1)
               XWRK2 = XHIGT (2)
               XWRK3 = XHIGT (IFIN)
               if (XWRK2 < XWRK1) then
                  XHIGT (1) = XWRK2
                  XHIGT (2) = XWRK1
                  XWRK2 = XWRK1
               end if
               if (XWRK2 > XWRK3) then
                  XHIGT (IFIN) = XWRK2
                  XHIGT (2) = XWRK3
                  XWRK2 = XWRK3
                  if (XWRK2 < XHIGT(1)) then
                     XHIGT (2) = XHIGT (1)
                     XHIGT (1) = XWRK2
                  end if
               end if
!
               XWRK1 = XHIGT (1)
               JLOW = JLOW + 1
               XLOWT (JLOW) = XWRK1
               XPIV = XWRK1 + 0.5 * (XHIGT(IFIN)-XWRK1)
!
!  One takes values <= pivot to XLOWT
!  Again, 2 parts, one where we take care of the remaining
!  high values because we might still need them, and the
!  other when we know that we will have more than enough
!  low values in the end.
!
               JHIG = 0
               do ICRS = 2, IFIN
                  if (XHIGT(ICRS) <= XPIV) then
                     JLOW = JLOW + 1
                     XLOWT (JLOW) = XHIGT (ICRS)
                     if (JLOW >= INTH) exit
                  else
                     JHIG = JHIG + 1
                     XHIGT (JHIG) = XHIGT (ICRS)
                  end if
               end do
!
               do ICRS = ICRS + 1, IFIN
                  if (XHIGT(ICRS) <= XPIV) then
                     JLOW = JLOW + 1
                     XLOWT (JLOW) = XHIGT (ICRS)
                  end if
               end do
            end select
!
!
         case (1)
!
!  Only 1 value is missing in low part
!
            XMIN = XHIGT(1)
            IHIG = 1
            do ICRS = 2, JHIG
               if (XHIGT(ICRS) < XMIN) then
                  XMIN = XHIGT(ICRS)
                  IHIG = ICRS
               end if
            end do
!
            VALNTH = XHIGT (IHIG)
            return
!
!
         case (0)
!
!  Low part is exactly what we want
!
            exit
!
!
         case (-5:-1)
!
!  Only few values too many in low part
!
            XHIGT (1) = XLOWT (1)
            ILOW = 1 + INTH - JLOW
            do ICRS = 2, INTH
               XWRK = XLOWT (ICRS)
               do IDCR = ICRS - 1, MAX (1, ILOW), - 1
                  if (XWRK < XHIGT(IDCR)) then
                     XHIGT (IDCR+1) = XHIGT (IDCR)
                  else
                     exit
                  end if
               end do
               XHIGT (IDCR+1) = XWRK
               ILOW = ILOW + 1
            end do
!
            XWRK1 = XHIGT(INTH)
            ILOW = 2*INTH - JLOW
            do ICRS = INTH + 1, JLOW
               if (XLOWT (ICRS) < XWRK1) then
                  XWRK = XLOWT (ICRS)
                  do IDCR = INTH - 1, MAX (1, ILOW), - 1
                     if (XWRK >= XHIGT(IDCR)) exit
                     XHIGT (IDCR+1) = XHIGT (IDCR)
                  end do
                  XHIGT (IDCR+1) = XLOWT (ICRS)
                  XWRK1 = XHIGT(INTH)
               end if
               ILOW = ILOW + 1
            end do
!
            VALNTH = XHIGT(INTH)
            return
!
!
         case (:-6)
!
! last case: too many values in low part
!

            IMIL = (JLOW+1) / 2
            IFIN = JLOW
!
!  One chooses a pivot from 1st, last, and middle values
!
            if (XLOWT(IMIL) < XLOWT(1)) then
               XWRK = XLOWT (1)
               XLOWT (1) = XLOWT (IMIL)
               XLOWT (IMIL) = XWRK
            end if
            if (XLOWT(IMIL) > XLOWT(IFIN)) then
               XWRK = XLOWT (IFIN)
               XLOWT (IFIN) = XLOWT (IMIL)
               XLOWT (IMIL) = XWRK
               if (XLOWT(IMIL) < XLOWT(1)) then
                  XWRK = XLOWT (1)
                  XLOWT (1) = XLOWT (IMIL)
                  XLOWT (IMIL) = XWRK
               end if
            end if
            if (IFIN <= 3) exit
!
            XPIV = XLOWT(1) + REAL(INTH)/REAL(JLOW+INTH) * &
                              (XLOWT(IFIN)-XLOWT(1))

!
!  One takes values > XPIV to XHIGT
!
            JHIG = 0
            JLOW = 0
!
            if (XLOWT(IFIN) > XPIV) then
               ICRS = 0
               do
                  ICRS = ICRS + 1
                  if (XLOWT(ICRS) > XPIV) then
                     JHIG = JHIG + 1
                     XHIGT (JHIG) = XLOWT (ICRS)
                     if (ICRS >= IFIN) exit
                  else
                     JLOW = JLOW + 1
                     XLOWT (JLOW) = XLOWT (ICRS)
                     if (JLOW >= INTH) exit
                  end if
               end do
!
               if (ICRS < IFIN) then
                  do
                     ICRS = ICRS + 1
                     if (XLOWT(ICRS) <= XPIV) then
                        JLOW = JLOW + 1
                        XLOWT (JLOW) = XLOWT (ICRS)
                     else
                        if (ICRS >= IFIN) exit
                     end if
                  end do
               end if
            else
               do ICRS = 1, IFIN
                  if (XLOWT(ICRS) > XPIV) then
                     JHIG = JHIG + 1
                     XHIGT (JHIG) = XLOWT (ICRS)
                  else
                     JLOW = JLOW + 1
                     XLOWT (JLOW) = XLOWT (ICRS)
                     if (JLOW >= INTH) exit
                  end if
               end do
!
               do ICRS = ICRS + 1, IFIN
                  if (XLOWT(ICRS) <= XPIV) then
                     JLOW = JLOW + 1
                     XLOWT (JLOW) = XLOWT (ICRS)
                  end if
               end do
            end if
!
         end select
!
      end do
!
!  Now, we only need to find maximum of the 1:INTH set
!
      VALNTH = MAXVAL (XLOWT (1:INTH))
      return
!
!
end function D_valnth

pure function R_valnth (XDONT, NORD) result (valnth)
!  Return NORDth value of XDONT, i.e fractile of order NORD/SIZE(XDONT).
! __________________________________________________________
!  This routine uses a pivoting strategy such as the one of
!  finding the median based on the quicksort algorithm, but
!  we skew the pivot choice to try to bring it to NORD as
!  fast as possible. It uses 2 temporary arrays, where it
!  stores the indices of the values smaller than the pivot
!  (ILOWT), and the indices of values larger than the pivot
!  that we might still need later on (IHIGT). It iterates
!  until it can bring the number of values in ILOWT to
!  exactly NORD, and then finds the maximum of this set.
!  Michel Olagnon - Aug. 2000
! __________________________________________________________
! _________________________________________________________
      real, dimension (:), intent (in) :: XDONT
      real :: valnth
      integer, intent (in) :: NORD
! __________________________________________________________
      real, dimension (SIZE(XDONT)) :: XLOWT, XHIGT
      real :: XPIV, XPIV0, XWRK, XWRK1, XWRK2, XWRK3, XMIN, XMAX
!
      integer :: NDON, JHIG, JLOW, IHIG
      integer :: IMIL, IFIN, ICRS, IDCR, ILOW
      integer :: JLM2, JLM1, JHM2, JHM1, INTH
!
      NDON = SIZE (XDONT)
      INTH = MAX (MIN (NORD, NDON), 1)
!
!    First loop is used to fill-in XLOWT, XHIGT at the same time
!
      if (NDON < 2) then
         if (INTH == 1) VALNTH = XDONT (1)
         return
      end if
!
!  One chooses a pivot, best estimate possible to put fractile near
!  mid-point of the set of low values.
!
      if (XDONT(2) < XDONT(1)) then
         XLOWT (1) = XDONT(2)
         XHIGT (1) = XDONT(1)
      else
         XLOWT (1) = XDONT(1)
         XHIGT (1) = XDONT(2)
      end if
!
      if (NDON < 3) then
         if (INTH == 1) VALNTH = XLOWT (1)
         if (INTH == 2) VALNTH = XHIGT (1)
         return
      end if
!
      if (XDONT(3) < XHIGT(1)) then
         XHIGT (2) = XHIGT (1)
         if (XDONT(3) < XLOWT(1)) then
            XHIGT (1) = XLOWT (1)
            XLOWT (1) = XDONT(3)
         else
            XHIGT (1) = XDONT(3)
         end if
      else
         XHIGT (2) = XDONT(3)
      end if
!
      if (NDON < 4) then
         if (INTH == 1) then
             VALNTH = XLOWT (1)
         else
             VALNTH = XHIGT (INTH - 1)
         end if
         return
      end if
!
      if (XDONT(NDON) < XHIGT(1)) then
         XHIGT (3) = XHIGT (2)
         XHIGT (2) = XHIGT (1)
         if (XDONT(NDON) < XLOWT(1)) then
            XHIGT (1) = XLOWT (1)
            XLOWT (1) = XDONT(NDON)
         else
            XHIGT (1) = XDONT(NDON)
         end if
      else
         XHIGT (3) = XDONT(NDON)
      end if
!
      if (NDON < 5) then
         if (INTH == 1) then
             VALNTH = XLOWT (1)
         else
             VALNTH = XHIGT (INTH - 1)
         end if
         return
      end if
!

      JLOW = 1
      JHIG = 3
      XPIV = XLOWT(1) + REAL(2*INTH)/REAL(NDON+INTH) * (XHIGT(3)-XLOWT(1))
      if (XPIV >= XHIGT(1)) then
         XPIV = XLOWT(1) + REAL(2*INTH)/REAL(NDON+INTH) * &
                           (XHIGT(2)-XLOWT(1))
         if (XPIV >= XHIGT(1)) &
             XPIV = XLOWT(1) + REAL (2*INTH) / REAL (NDON+INTH) * &
                               (XHIGT(1)-XLOWT(1))
      end if
      XPIV0 = XPIV
!
!  One puts values > pivot in the end and those <= pivot
!  at the beginning. This is split in 2 cases, so that
!  we can skip the loop test a number of times.
!  As we are also filling in the work arrays at the same time
!  we stop filling in the XHIGT array as soon as we have more
!  than enough values in XLOWT.
!
!
      if (XDONT(NDON) > XPIV) then
         ICRS = 3
         do
            ICRS = ICRS + 1
            if (XDONT(ICRS) > XPIV) then
               if (ICRS >= NDON) exit
               JHIG = JHIG + 1
               XHIGT (JHIG) = XDONT(ICRS)
            else
               JLOW = JLOW + 1
               XLOWT (JLOW) = XDONT(ICRS)
               if (JLOW >= INTH) exit
            end if
         end do
!
!  One restricts further processing because it is no use
!  to store more high values
!
         if (ICRS < NDON-1) then
            do
               ICRS = ICRS + 1
               if (XDONT(ICRS) <= XPIV) then
                  JLOW = JLOW + 1
                  XLOWT (JLOW) = XDONT(ICRS)
               else if (ICRS >= NDON) then
                  exit
               end if
            end do
         end if
!
!
      else
!
!  Same as above, but this is not as easy to optimize, so the
!  DO-loop is kept
!
         do ICRS = 4, NDON - 1
            if (XDONT(ICRS) > XPIV) then
               JHIG = JHIG + 1
               XHIGT (JHIG) = XDONT(ICRS)
            else
               JLOW = JLOW + 1
               XLOWT (JLOW) = XDONT(ICRS)
               if (JLOW >= INTH) exit
            end if
         end do
!
         if (ICRS < NDON-1) then
            do
               ICRS = ICRS + 1
               if (XDONT(ICRS) <= XPIV) then
                  if (ICRS >= NDON) exit
                  JLOW = JLOW + 1
                  XLOWT (JLOW) = XDONT(ICRS)
               end if
            end do
         end if
      end if
!
      JLM2 = 0
      JLM1 = 0
      JHM2 = 0
      JHM1 = 0
      do
         if (JLM2 == JLOW .and. JHM2 == JHIG) then
!
!   We are oscillating. Perturbate by bringing JLOW closer by one
!   to INTH
!
             if (INTH > JLOW) then
                XMIN = XHIGT(1)
                IHIG = 1
                do ICRS = 2, JHIG
                   if (XHIGT(ICRS) < XMIN) then
                      XMIN = XHIGT(ICRS)
                      IHIG = ICRS
                   end if
                end do
!
                JLOW = JLOW + 1
                XLOWT (JLOW) = XHIGT (IHIG)
                XHIGT (IHIG) = XHIGT (JHIG)
                JHIG = JHIG - 1
             else

                XMAX = XLOWT (JLOW)
                JLOW = JLOW - 1
                do ICRS = 1, JLOW
                   if (XLOWT(ICRS) > XMAX) then
                      XWRK = XMAX
                      XMAX = XLOWT(ICRS)
                      XLOWT (ICRS) = XWRK
                   end if
                end do
             end if
         end if
         JLM2 = JLM1
         JLM1 = JLOW
         JHM2 = JHM1
         JHM1 = JHIG
!
!   We try to bring the number of values in the low values set
!   closer to INTH.
!
         select case (INTH-JLOW)
         case (2:)
!
!   Not enough values in low part, at least 2 are missing
!
            INTH = INTH - JLOW
            JLOW = 0
            select case (JHIG)
!!!!!           CASE DEFAULT
!!!!!              write (unit=*,fmt=*) "Assertion failed"
!!!!!              STOP
!
!   We make a special case when we have so few values in
!   the high values set that it is bad performance to choose a pivot
!   and apply the general algorithm.
!
            case (2)
               if (XHIGT(1) <= XHIGT(2)) then
                  JLOW = JLOW + 1
                  XLOWT (JLOW) = XHIGT (1)
                  JLOW = JLOW + 1
                  XLOWT (JLOW) = XHIGT (2)
               else
                  JLOW = JLOW + 1
                  XLOWT (JLOW) = XHIGT (2)
                  JLOW = JLOW + 1
                  XLOWT (JLOW) = XHIGT (1)
               end if
               exit
!
            case (3)
!
!
               XWRK1 = XHIGT (1)
               XWRK2 = XHIGT (2)
               XWRK3 = XHIGT (3)
               if (XWRK2 < XWRK1) then
                  XHIGT (1) = XWRK2
                  XHIGT (2) = XWRK1
                  XWRK2 = XWRK1
               end if
               if (XWRK2 > XWRK3) then
                  XHIGT (3) = XWRK2
                  XHIGT (2) = XWRK3
                  XWRK2 = XWRK3
                  if (XWRK2 < XHIGT(1)) then
                     XHIGT (2) = XHIGT (1)
                     XHIGT (1) = XWRK2
                  end if
               end if
               JHIG = 0
               do ICRS = JLOW + 1, INTH
                  JHIG = JHIG + 1
                  XLOWT (ICRS) = XHIGT (JHIG)
               end do
               JLOW = INTH
               exit
!
            case (4:)
!
!
               XPIV0 = XPIV
               IFIN = JHIG
!
!  One chooses a pivot from the 2 first values and the last one.
!  This should ensure sufficient renewal between iterations to
!  avoid worst case behavior effects.
!
               XWRK1 = XHIGT (1)
               XWRK2 = XHIGT (2)
               XWRK3 = XHIGT (IFIN)
               if (XWRK2 < XWRK1) then
                  XHIGT (1) = XWRK2
                  XHIGT (2) = XWRK1
                  XWRK2 = XWRK1
               end if
               if (XWRK2 > XWRK3) then
                  XHIGT (IFIN) = XWRK2
                  XHIGT (2) = XWRK3
                  XWRK2 = XWRK3
                  if (XWRK2 < XHIGT(1)) then
                     XHIGT (2) = XHIGT (1)
                     XHIGT (1) = XWRK2
                  end if
               end if
!
               XWRK1 = XHIGT (1)
               JLOW = JLOW + 1
               XLOWT (JLOW) = XWRK1
               XPIV = XWRK1 + 0.5 * (XHIGT(IFIN)-XWRK1)
!
!  One takes values <= pivot to XLOWT
!  Again, 2 parts, one where we take care of the remaining
!  high values because we might still need them, and the
!  other when we know that we will have more than enough
!  low values in the end.
!
               JHIG = 0
               do ICRS = 2, IFIN
                  if (XHIGT(ICRS) <= XPIV) then
                     JLOW = JLOW + 1
                     XLOWT (JLOW) = XHIGT (ICRS)
                     if (JLOW >= INTH) exit
                  else
                     JHIG = JHIG + 1
                     XHIGT (JHIG) = XHIGT (ICRS)
                  end if
               end do
!
               do ICRS = ICRS + 1, IFIN
                  if (XHIGT(ICRS) <= XPIV) then
                     JLOW = JLOW + 1
                     XLOWT (JLOW) = XHIGT (ICRS)
                  end if
               end do
            end select
!
!
         case (1)
!
!  Only 1 value is missing in low part
!
            XMIN = XHIGT(1)
            IHIG = 1
            do ICRS = 2, JHIG
               if (XHIGT(ICRS) < XMIN) then
                  XMIN = XHIGT(ICRS)
                  IHIG = ICRS
               end if
            end do
!
            VALNTH = XHIGT (IHIG)
            return
!
!
         case (0)
!
!  Low part is exactly what we want
!
            exit
!
!
         case (-5:-1)
!
!  Only few values too many in low part
!
            XHIGT (1) = XLOWT (1)
            ILOW = 1 + INTH - JLOW
            do ICRS = 2, INTH
               XWRK = XLOWT (ICRS)
               do IDCR = ICRS - 1, MAX (1, ILOW), - 1
                  if (XWRK < XHIGT(IDCR)) then
                     XHIGT (IDCR+1) = XHIGT (IDCR)
                  else
                     exit
                  end if
               end do
               XHIGT (IDCR+1) = XWRK
               ILOW = ILOW + 1
            end do
!
            XWRK1 = XHIGT(INTH)
            ILOW = 2*INTH - JLOW
            do ICRS = INTH + 1, JLOW
               if (XLOWT (ICRS) < XWRK1) then
                  XWRK = XLOWT (ICRS)
                  do IDCR = INTH - 1, MAX (1, ILOW), - 1
                     if (XWRK >= XHIGT(IDCR)) exit
                     XHIGT (IDCR+1) = XHIGT (IDCR)
                  end do
                  XHIGT (IDCR+1) = XLOWT (ICRS)
                  XWRK1 = XHIGT(INTH)
               end if
               ILOW = ILOW + 1
            end do
!
            VALNTH = XHIGT(INTH)
            return
!
!
         case (:-6)
!
! last case: too many values in low part
!

            IMIL = (JLOW+1) / 2
            IFIN = JLOW
!
!  One chooses a pivot from 1st, last, and middle values
!
            if (XLOWT(IMIL) < XLOWT(1)) then
               XWRK = XLOWT (1)
               XLOWT (1) = XLOWT (IMIL)
               XLOWT (IMIL) = XWRK
            end if
            if (XLOWT(IMIL) > XLOWT(IFIN)) then
               XWRK = XLOWT (IFIN)
               XLOWT (IFIN) = XLOWT (IMIL)
               XLOWT (IMIL) = XWRK
               if (XLOWT(IMIL) < XLOWT(1)) then
                  XWRK = XLOWT (1)
                  XLOWT (1) = XLOWT (IMIL)
                  XLOWT (IMIL) = XWRK
               end if
            end if
            if (IFIN <= 3) exit
!
            XPIV = XLOWT(1) + REAL(INTH)/REAL(JLOW+INTH) * &
                              (XLOWT(IFIN)-XLOWT(1))

!
!  One takes values > XPIV to XHIGT
!
            JHIG = 0
            JLOW = 0
!
            if (XLOWT(IFIN) > XPIV) then
               ICRS = 0
               do
                  ICRS = ICRS + 1
                  if (XLOWT(ICRS) > XPIV) then
                     JHIG = JHIG + 1
                     XHIGT (JHIG) = XLOWT (ICRS)
                     if (ICRS >= IFIN) exit
                  else
                     JLOW = JLOW + 1
                     XLOWT (JLOW) = XLOWT (ICRS)
                     if (JLOW >= INTH) exit
                  end if
               end do
!
               if (ICRS < IFIN) then
                  do
                     ICRS = ICRS + 1
                     if (XLOWT(ICRS) <= XPIV) then
                        JLOW = JLOW + 1
                        XLOWT (JLOW) = XLOWT (ICRS)
                     else
                        if (ICRS >= IFIN) exit
                     end if
                  end do
               end if
            else
               do ICRS = 1, IFIN
                  if (XLOWT(ICRS) > XPIV) then
                     JHIG = JHIG + 1
                     XHIGT (JHIG) = XLOWT (ICRS)
                  else
                     JLOW = JLOW + 1
                     XLOWT (JLOW) = XLOWT (ICRS)
                     if (JLOW >= INTH) exit
                  end if
               end do
!
               do ICRS = ICRS + 1, IFIN
                  if (XLOWT(ICRS) <= XPIV) then
                     JLOW = JLOW + 1
                     XLOWT (JLOW) = XLOWT (ICRS)
                  end if
               end do
            end if
!
         end select
!
      end do
!
!  Now, we only need to find maximum of the 1:INTH set
!
      VALNTH = MAXVAL (XLOWT (1:INTH))
      return
!
!
end function R_valnth
pure function I_valnth (XDONT, NORD) result (valnth)
!  Return NORDth value of XDONT, i.e fractile of order NORD/SIZE(XDONT).
! __________________________________________________________
!  This routine uses a pivoting strategy such as the one of
!  finding the median based on the quicksort algorithm, but
!  we skew the pivot choice to try to bring it to NORD as
!  fast as possible. It uses 2 temporary arrays, where it
!  stores the indices of the values smaller than the pivot
!  (ILOWT), and the indices of values larger than the pivot
!  that we might still need later on (IHIGT). It iterates
!  until it can bring the number of values in ILOWT to
!  exactly NORD, and then finds the maximum of this set.
!  Michel Olagnon - Aug. 2000
! __________________________________________________________
! __________________________________________________________
      integer, dimension (:), intent (in) :: XDONT
      integer :: valnth
      integer, intent (in) :: NORD
! __________________________________________________________
      integer, dimension (SIZE(XDONT)) :: XLOWT, XHIGT
      integer :: XPIV, XPIV0, XWRK, XWRK1, XWRK2, XWRK3, XMIN, XMAX
!
      integer :: NDON, JHIG, JLOW, IHIG
      integer :: IMIL, IFIN, ICRS, IDCR, ILOW
      integer :: JLM2, JLM1, JHM2, JHM1, INTH
!
      NDON = SIZE (XDONT)
      INTH = MAX (MIN (NORD, NDON), 1)
!
!    First loop is used to fill-in XLOWT, XHIGT at the same time
!
      if (NDON < 2) then
         if (INTH == 1) VALNTH = XDONT (1)
         return
      end if
!
!  One chooses a pivot, best estimate possible to put fractile near
!  mid-point of the set of low values.
!
      if (XDONT(2) < XDONT(1)) then
         XLOWT (1) = XDONT(2)
         XHIGT (1) = XDONT(1)
      else
         XLOWT (1) = XDONT(1)
         XHIGT (1) = XDONT(2)
      end if
!
      if (NDON < 3) then
         if (INTH == 1) VALNTH = XLOWT (1)
         if (INTH == 2) VALNTH = XHIGT (1)
         return
      end if
!
      if (XDONT(3) < XHIGT(1)) then
         XHIGT (2) = XHIGT (1)
         if (XDONT(3) < XLOWT(1)) then
            XHIGT (1) = XLOWT (1)
            XLOWT (1) = XDONT(3)
         else
            XHIGT (1) = XDONT(3)
         end if
      else
         XHIGT (2) = XDONT(3)
      end if
!
      if (NDON < 4) then
         if (INTH == 1) then
             VALNTH = XLOWT (1)
         else
             VALNTH = XHIGT (INTH - 1)
         end if
         return
      end if
!
      if (XDONT(NDON) < XHIGT(1)) then
         XHIGT (3) = XHIGT (2)
         XHIGT (2) = XHIGT (1)
         if (XDONT(NDON) < XLOWT(1)) then
            XHIGT (1) = XLOWT (1)
            XLOWT (1) = XDONT(NDON)
         else
            XHIGT (1) = XDONT(NDON)
         end if
      else
         XHIGT (3) = XDONT(NDON)
      end if
!
      if (NDON < 5) then
         if (INTH == 1) then
             VALNTH = XLOWT (1)
         else
             VALNTH = XHIGT (INTH - 1)
         end if
         return
      end if
!

      JLOW = 1
      JHIG = 3
      XPIV = XLOWT(1) + REAL(2*INTH)/REAL(NDON+INTH) * (XHIGT(3)-XLOWT(1))
      if (XPIV >= XHIGT(1)) then
         XPIV = XLOWT(1) + REAL(2*INTH)/REAL(NDON+INTH) * &
                           (XHIGT(2)-XLOWT(1))
         if (XPIV >= XHIGT(1)) &
             XPIV = XLOWT(1) + REAL (2*INTH) / REAL (NDON+INTH) * &
                               (XHIGT(1)-XLOWT(1))
      end if
      XPIV0 = XPIV
!
!  One puts values > pivot in the end and those <= pivot
!  at the beginning. This is split in 2 cases, so that
!  we can skip the loop test a number of times.
!  As we are also filling in the work arrays at the same time
!  we stop filling in the XHIGT array as soon as we have more
!  than enough values in XLOWT.
!
!
      if (XDONT(NDON) > XPIV) then
         ICRS = 3
         do
            ICRS = ICRS + 1
            if (XDONT(ICRS) > XPIV) then
               if (ICRS >= NDON) exit
               JHIG = JHIG + 1
               XHIGT (JHIG) = XDONT(ICRS)
            else
               JLOW = JLOW + 1
               XLOWT (JLOW) = XDONT(ICRS)
               if (JLOW >= INTH) exit
            end if
         end do
!
!  One restricts further processing because it is no use
!  to store more high values
!
         if (ICRS < NDON-1) then
            do
               ICRS = ICRS + 1
               if (XDONT(ICRS) <= XPIV) then
                  JLOW = JLOW + 1
                  XLOWT (JLOW) = XDONT(ICRS)
               else if (ICRS >= NDON) then
                  exit
               end if
            end do
         end if
!
!
      else
!
!  Same as above, but this is not as easy to optimize, so the
!  DO-loop is kept
!
         do ICRS = 4, NDON - 1
            if (XDONT(ICRS) > XPIV) then
               JHIG = JHIG + 1
               XHIGT (JHIG) = XDONT(ICRS)
            else
               JLOW = JLOW + 1
               XLOWT (JLOW) = XDONT(ICRS)
               if (JLOW >= INTH) exit
            end if
         end do
!
         if (ICRS < NDON-1) then
            do
               ICRS = ICRS + 1
               if (XDONT(ICRS) <= XPIV) then
                  if (ICRS >= NDON) exit
                  JLOW = JLOW + 1
                  XLOWT (JLOW) = XDONT(ICRS)
               end if
            end do
         end if
      end if
!
      JLM2 = 0
      JLM1 = 0
      JHM2 = 0
      JHM1 = 0
      do
         if (JLM2 == JLOW .and. JHM2 == JHIG) then
!
!   We are oscillating. Perturbate by bringing JLOW closer by one
!   to INTH
!
             if (INTH > JLOW) then
                XMIN = XHIGT(1)
                IHIG = 1
                do ICRS = 2, JHIG
                   if (XHIGT(ICRS) < XMIN) then
                      XMIN = XHIGT(ICRS)
                      IHIG = ICRS
                   end if
                end do
!
                JLOW = JLOW + 1
                XLOWT (JLOW) = XHIGT (IHIG)
                XHIGT (IHIG) = XHIGT (JHIG)
                JHIG = JHIG - 1
             else

                XMAX = XLOWT (JLOW)
                JLOW = JLOW - 1
                do ICRS = 1, JLOW
                   if (XLOWT(ICRS) > XMAX) then
                      XWRK = XMAX
                      XMAX = XLOWT(ICRS)
                      XLOWT (ICRS) = XWRK
                   end if
                end do
             end if
         end if
         JLM2 = JLM1
         JLM1 = JLOW
         JHM2 = JHM1
         JHM1 = JHIG
!
!   We try to bring the number of values in the low values set
!   closer to INTH.
!
         select case (INTH-JLOW)
         case (2:)
!
!   Not enough values in low part, at least 2 are missing
!
            INTH = INTH - JLOW
            JLOW = 0
            select case (JHIG)
!!!!!           CASE DEFAULT
!!!!!              write (unit=*,fmt=*) "Assertion failed"
!!!!!              STOP
!
!   We make a special case when we have so few values in
!   the high values set that it is bad performance to choose a pivot
!   and apply the general algorithm.
!
            case (2)
               if (XHIGT(1) <= XHIGT(2)) then
                  JLOW = JLOW + 1
                  XLOWT (JLOW) = XHIGT (1)
                  JLOW = JLOW + 1
                  XLOWT (JLOW) = XHIGT (2)
               else
                  JLOW = JLOW + 1
                  XLOWT (JLOW) = XHIGT (2)
                  JLOW = JLOW + 1
                  XLOWT (JLOW) = XHIGT (1)
               end if
               exit
!
            case (3)
!
!
               XWRK1 = XHIGT (1)
               XWRK2 = XHIGT (2)
               XWRK3 = XHIGT (3)
               if (XWRK2 < XWRK1) then
                  XHIGT (1) = XWRK2
                  XHIGT (2) = XWRK1
                  XWRK2 = XWRK1
               end if
               if (XWRK2 > XWRK3) then
                  XHIGT (3) = XWRK2
                  XHIGT (2) = XWRK3
                  XWRK2 = XWRK3
                  if (XWRK2 < XHIGT(1)) then
                     XHIGT (2) = XHIGT (1)
                     XHIGT (1) = XWRK2
                  end if
               end if
               JHIG = 0
               do ICRS = JLOW + 1, INTH
                  JHIG = JHIG + 1
                  XLOWT (ICRS) = XHIGT (JHIG)
               end do
               JLOW = INTH
               exit
!
            case (4:)
!
!
               XPIV0 = XPIV
               IFIN = JHIG
!
!  One chooses a pivot from the 2 first values and the last one.
!  This should ensure sufficient renewal between iterations to
!  avoid worst case behavior effects.
!
               XWRK1 = XHIGT (1)
               XWRK2 = XHIGT (2)
               XWRK3 = XHIGT (IFIN)
               if (XWRK2 < XWRK1) then
                  XHIGT (1) = XWRK2
                  XHIGT (2) = XWRK1
                  XWRK2 = XWRK1
               end if
               if (XWRK2 > XWRK3) then
                  XHIGT (IFIN) = XWRK2
                  XHIGT (2) = XWRK3
                  XWRK2 = XWRK3
                  if (XWRK2 < XHIGT(1)) then
                     XHIGT (2) = XHIGT (1)
                     XHIGT (1) = XWRK2
                  end if
               end if
!
               XWRK1 = XHIGT (1)
               JLOW = JLOW + 1
               XLOWT (JLOW) = XWRK1
               XPIV = XWRK1 + 0.5 * (XHIGT(IFIN)-XWRK1)
!
!  One takes values <= pivot to XLOWT
!  Again, 2 parts, one where we take care of the remaining
!  high values because we might still need them, and the
!  other when we know that we will have more than enough
!  low values in the end.
!
               JHIG = 0
               do ICRS = 2, IFIN
                  if (XHIGT(ICRS) <= XPIV) then
                     JLOW = JLOW + 1
                     XLOWT (JLOW) = XHIGT (ICRS)
                     if (JLOW >= INTH) exit
                  else
                     JHIG = JHIG + 1
                     XHIGT (JHIG) = XHIGT (ICRS)
                  end if
               end do
!
               do ICRS = ICRS + 1, IFIN
                  if (XHIGT(ICRS) <= XPIV) then
                     JLOW = JLOW + 1
                     XLOWT (JLOW) = XHIGT (ICRS)
                  end if
               end do
            end select
!
!
         case (1)
!
!  Only 1 value is missing in low part
!
            XMIN = XHIGT(1)
            IHIG = 1
            do ICRS = 2, JHIG
               if (XHIGT(ICRS) < XMIN) then
                  XMIN = XHIGT(ICRS)
                  IHIG = ICRS
               end if
            end do
!
            VALNTH = XHIGT (IHIG)
            return
!
!
         case (0)
!
!  Low part is exactly what we want
!
            exit
!
!
         case (-5:-1)
!
!  Only few values too many in low part
!
            XHIGT (1) = XLOWT (1)
            ILOW = 1 + INTH - JLOW
            do ICRS = 2, INTH
               XWRK = XLOWT (ICRS)
               do IDCR = ICRS - 1, MAX (1, ILOW), - 1
                  if (XWRK < XHIGT(IDCR)) then
                     XHIGT (IDCR+1) = XHIGT (IDCR)
                  else
                     exit
                  end if
               end do
               XHIGT (IDCR+1) = XWRK
               ILOW = ILOW + 1
            end do
!
            XWRK1 = XHIGT(INTH)
            ILOW = 2*INTH - JLOW
            do ICRS = INTH + 1, JLOW
               if (XLOWT (ICRS) < XWRK1) then
                  XWRK = XLOWT (ICRS)
                  do IDCR = INTH - 1, MAX (1, ILOW), - 1
                     if (XWRK >= XHIGT(IDCR)) exit
                     XHIGT (IDCR+1) = XHIGT (IDCR)
                  end do
                  XHIGT (IDCR+1) = XLOWT (ICRS)
                  XWRK1 = XHIGT(INTH)
               end if
               ILOW = ILOW + 1
            end do
!
            VALNTH = XHIGT(INTH)
            return
!
!
         case (:-6)
!
! last case: too many values in low part
!

            IMIL = (JLOW+1) / 2
            IFIN = JLOW
!
!  One chooses a pivot from 1st, last, and middle values
!
            if (XLOWT(IMIL) < XLOWT(1)) then
               XWRK = XLOWT (1)
               XLOWT (1) = XLOWT (IMIL)
               XLOWT (IMIL) = XWRK
            end if
            if (XLOWT(IMIL) > XLOWT(IFIN)) then
               XWRK = XLOWT (IFIN)
               XLOWT (IFIN) = XLOWT (IMIL)
               XLOWT (IMIL) = XWRK
               if (XLOWT(IMIL) < XLOWT(1)) then
                  XWRK = XLOWT (1)
                  XLOWT (1) = XLOWT (IMIL)
                  XLOWT (IMIL) = XWRK
               end if
            end if
            if (IFIN <= 3) exit
!
            XPIV = XLOWT(1) + REAL(INTH)/REAL(JLOW+INTH) * &
                              (XLOWT(IFIN)-XLOWT(1))

!
!  One takes values > XPIV to XHIGT
!
            JHIG = 0
            JLOW = 0
!
            if (XLOWT(IFIN) > XPIV) then
               ICRS = 0
               do
                  ICRS = ICRS + 1
                  if (XLOWT(ICRS) > XPIV) then
                     JHIG = JHIG + 1
                     XHIGT (JHIG) = XLOWT (ICRS)
                     if (ICRS >= IFIN) exit
                  else
                     JLOW = JLOW + 1
                     XLOWT (JLOW) = XLOWT (ICRS)
                     if (JLOW >= INTH) exit
                  end if
               end do
!
               if (ICRS < IFIN) then
                  do
                     ICRS = ICRS + 1
                     if (XLOWT(ICRS) <= XPIV) then
                        JLOW = JLOW + 1
                        XLOWT (JLOW) = XLOWT (ICRS)
                     else
                        if (ICRS >= IFIN) exit
                     end if
                  end do
               end if
            else
               do ICRS = 1, IFIN
                  if (XLOWT(ICRS) > XPIV) then
                     JHIG = JHIG + 1
                     XHIGT (JHIG) = XLOWT (ICRS)
                  else
                     JLOW = JLOW + 1
                     XLOWT (JLOW) = XLOWT (ICRS)
                     if (JLOW >= INTH) exit
                  end if
               end do
!
               do ICRS = ICRS + 1, IFIN
                  if (XLOWT(ICRS) <= XPIV) then
                     JLOW = JLOW + 1
                     XLOWT (JLOW) = XLOWT (ICRS)
                  end if
               end do
            end if
!
         end select
!
      end do
!
!  Now, we only need to find maximum of the 1:INTH set
!
      VALNTH = MAXVAL (XLOWT (1:INTH))
      return
!
!
end function I_valnth

!-------------------------------------------------------------------------------
!*******************************************************************************
!*******************************************************************************
!-------------------------------------------------------------------------------


pure function D_indnth (XDONT, NORD) result (INDNTH)
!  Return NORDth value of XDONT, i.e fractile of order NORD/SIZE(XDONT).
! __________________________________________________________
!  This routine uses a pivoting strategy such as the one of
!  finding the median based on the quicksort algorithm, but
!  we skew the pivot choice to try to bring it to NORD as
!  fast as possible. It uses 2 temporary arrays, where it
!  stores the indices of the values smaller than the pivot
!  (ILOWT), and the indices of values larger than the pivot
!  that we might still need later on (IHIGT). It iterates
!  until it can bring the number of values in ILOWT to
!  exactly NORD, and then finds the maximum of this set.
!  Michel Olagnon - Aug. 2000
! __________________________________________________________
! __________________________________________________________
      real (kind=kdp), dimension (:), intent (in) :: XDONT
      integer :: INDNTH
      integer, intent (in) :: NORD
! __________________________________________________________
      real (kind=kdp) :: XPIV, XPIV0, XWRK, XWRK1, XMIN, XMAX
!
      integer, dimension (NORD) :: IRNGT
      integer, dimension (SIZE(XDONT)) :: ILOWT, IHIGT
      integer :: NDON, JHIG, JLOW, IHIG, IWRK, IWRK1, IWRK2, IWRK3
      integer :: IMIL, IFIN, ICRS, IDCR, ILOW
      integer :: JLM2, JLM1, JHM2, JHM1, INTH
!
      NDON = SIZE (XDONT)
      INTH = NORD
!
!    First loop is used to fill-in ILOWT, IHIGT at the same time
!
      if (NDON < 2) then
         if (INTH == 1) INDNTH = 1
         return
      end if
!
!  One chooses a pivot, best estimate possible to put fractile near
!  mid-point of the set of low values.
!
      if (XDONT(2) < XDONT(1)) then
         ILOWT (1) = 2
         IHIGT (1) = 1
      else
         ILOWT (1) = 1
         IHIGT (1) = 2
      end if
!
      if (NDON < 3) then
         if (INTH == 1) INDNTH = ILOWT (1)
         if (INTH == 2) INDNTH = IHIGT (1)
         return
      end if
!
      if (XDONT(3) < XDONT(IHIGT(1))) then
         IHIGT (2) = IHIGT (1)
         if (XDONT(3) < XDONT(ILOWT(1))) then
            IHIGT (1) = ILOWT (1)
            ILOWT (1) = 3
         else
            IHIGT (1) = 3
         end if
      else
         IHIGT (2) = 3
      end if
!
      if (NDON < 4) then
         if (INTH == 1) INDNTH = ILOWT (1)
         if (INTH == 2) INDNTH = IHIGT (1)
         if (INTH == 3) INDNTH = IHIGT (2)
         return
      end if
!
      if (XDONT(NDON) < XDONT(IHIGT(1))) then
         IHIGT (3) = IHIGT (2)
         IHIGT (2) = IHIGT (1)
         if (XDONT(NDON) < XDONT(ILOWT(1))) then
            IHIGT (1) = ILOWT (1)
            ILOWT (1) = NDON
         else
            IHIGT (1) = NDON
         end if
      else
         IHIGT (3) = NDON
      end if
!
      if (NDON < 5) then
         if (INTH == 1) INDNTH = ILOWT (1)
         if (INTH == 2) INDNTH = IHIGT (1)
         if (INTH == 3) INDNTH = IHIGT (2)
         if (INTH == 4) INDNTH = IHIGT (3)
         return
      end if
!

      JLOW = 1
      JHIG = 3
      XPIV = XDONT (ILOWT(1)) + REAL(2*INTH)/REAL(NDON+INTH) * &
                                   (XDONT(IHIGT(3))-XDONT(ILOWT(1)))
      if (XPIV >= XDONT(IHIGT(1))) then
         XPIV = XDONT (ILOWT(1)) + REAL(2*INTH)/REAL(NDON+INTH) * &
                                      (XDONT(IHIGT(2))-XDONT(ILOWT(1)))
         if (XPIV >= XDONT(IHIGT(1))) &
             XPIV = XDONT (ILOWT(1)) + REAL (2*INTH) / REAL (NDON+INTH) * &
                                          (XDONT(IHIGT(1))-XDONT(ILOWT(1)))
      end if
      XPIV0 = XPIV
!
!  One puts values > pivot in the end and those <= pivot
!  at the beginning. This is split in 2 cases, so that
!  we can skip the loop test a number of times.
!  As we are also filling in the work arrays at the same time
!  we stop filling in the IHIGT array as soon as we have more
!  than enough values in ILOWT.
!
!
      if (XDONT(NDON) > XPIV) then
         ICRS = 3
         do
            ICRS = ICRS + 1
            if (XDONT(ICRS) > XPIV) then
               if (ICRS >= NDON) exit
               JHIG = JHIG + 1
               IHIGT (JHIG) = ICRS
            else
               JLOW = JLOW + 1
               ILOWT (JLOW) = ICRS
               if (JLOW >= INTH) exit
            end if
         end do
!
!  One restricts further processing because it is no use
!  to store more high values
!
         if (ICRS < NDON-1) then
            do
               ICRS = ICRS + 1
               if (XDONT(ICRS) <= XPIV) then
                  JLOW = JLOW + 1
                  ILOWT (JLOW) = ICRS
               else if (ICRS >= NDON) then
                  exit
               end if
            end do
         end if
!
!
      else
!
!  Same as above, but this is not as easy to optimize, so the
!  DO-loop is kept
!
         do ICRS = 4, NDON - 1
            if (XDONT(ICRS) > XPIV) then
               JHIG = JHIG + 1
               IHIGT (JHIG) = ICRS
            else
               JLOW = JLOW + 1
               ILOWT (JLOW) = ICRS
               if (JLOW >= INTH) exit
            end if
         end do
!
         if (ICRS < NDON-1) then
            do
               ICRS = ICRS + 1
               if (XDONT(ICRS) <= XPIV) then
                  if (ICRS >= NDON) exit
                  JLOW = JLOW + 1
                  ILOWT (JLOW) = ICRS
               end if
            end do
         end if
      end if
!
      JLM2 = 0
      JLM1 = 0
      JHM2 = 0
      JHM1 = 0
      do
         if (JLM2 == JLOW .and. JHM2 == JHIG) then
!
!   We are oscillating. Perturbate by bringing JLOW closer by one
!   to INTH
!
             if (INTH > JLOW) then
                XMIN = XDONT (IHIGT(1))
                IHIG = 1
                do ICRS = 2, JHIG
                   if (XDONT(IHIGT(ICRS)) < XMIN) then
                      XMIN = XDONT (IHIGT(ICRS))
                      IHIG = ICRS
                   end if
                end do
!
                JLOW = JLOW + 1
                ILOWT (JLOW) = IHIGT (IHIG)
                IHIGT (IHIG) = IHIGT (JHIG)
                JHIG = JHIG - 1
             else

                ILOW = ILOWT (1)
                XMAX = XDONT (ILOW)
                do ICRS = 2, JLOW
                   if (XDONT(ILOWT(ICRS)) > XMAX) then
                      IWRK = ILOWT (ICRS)
                      XMAX = XDONT (IWRK)
                      ILOWT (ICRS) = ILOW
                      ILOW = IWRK
                   end if
                end do
                JLOW = JLOW - 1
             end if
         end if
         JLM2 = JLM1
         JLM1 = JLOW
         JHM2 = JHM1
         JHM1 = JHIG
!
!   We try to bring the number of values in the low values set
!   closer to INTH.
!
         select case (INTH-JLOW)
         case (2:)
!
!   Not enough values in low part, at least 2 are missing
!
            INTH = INTH - JLOW
            JLOW = 0
            select case (JHIG)
!!!!!           CASE DEFAULT
!!!!!              write (unit=*,fmt=*) "Assertion failed"
!!!!!              STOP
!
!   We make a special case when we have so few values in
!   the high values set that it is bad performance to choose a pivot
!   and apply the general algorithm.
!
            case (2)
               if (XDONT(IHIGT(1)) <= XDONT(IHIGT(2))) then
                  JLOW = JLOW + 1
                  ILOWT (JLOW) = IHIGT (1)
                  JLOW = JLOW + 1
                  ILOWT (JLOW) = IHIGT (2)
               else
                  JLOW = JLOW + 1
                  ILOWT (JLOW) = IHIGT (2)
                  JLOW = JLOW + 1
                  ILOWT (JLOW) = IHIGT (1)
               end if
               exit
!
            case (3)
!
!
               IWRK1 = IHIGT (1)
               IWRK2 = IHIGT (2)
               IWRK3 = IHIGT (3)
               if (XDONT(IWRK2) < XDONT(IWRK1)) then
                  IHIGT (1) = IWRK2
                  IHIGT (2) = IWRK1
                  IWRK2 = IWRK1
               end if
               if (XDONT(IWRK2) > XDONT(IWRK3)) then
                  IHIGT (3) = IWRK2
                  IHIGT (2) = IWRK3
                  IWRK2 = IWRK3
                  if (XDONT(IWRK2) < XDONT(IHIGT(1))) then
                     IHIGT (2) = IHIGT (1)
                     IHIGT (1) = IWRK2
                  end if
               end if
               JHIG = 0
               do ICRS = JLOW + 1, INTH
                  JHIG = JHIG + 1
                  ILOWT (ICRS) = IHIGT (JHIG)
               end do
               JLOW = INTH
               exit
!
            case (4:)
!
!
               XPIV0 = XPIV
               IFIN = JHIG
!
!  One chooses a pivot from the 2 first values and the last one.
!  This should ensure sufficient renewal between iterations to
!  avoid worst case behavior effects.
!
               IWRK1 = IHIGT (1)
               IWRK2 = IHIGT (2)
               IWRK3 = IHIGT (IFIN)
               if (XDONT(IWRK2) < XDONT(IWRK1)) then
                  IHIGT (1) = IWRK2
                  IHIGT (2) = IWRK1
                  IWRK2 = IWRK1
               end if
               if (XDONT(IWRK2) > XDONT(IWRK3)) then
                  IHIGT (IFIN) = IWRK2
                  IHIGT (2) = IWRK3
                  IWRK2 = IWRK3
                  if (XDONT(IWRK2) < XDONT(IHIGT(1))) then
                     IHIGT (2) = IHIGT (1)
                     IHIGT (1) = IWRK2
                  end if
               end if
!
               IWRK1 = IHIGT (1)
               JLOW = JLOW + 1
               ILOWT (JLOW) = IWRK1
               XPIV = XDONT (IWRK1) + 0.5 * (XDONT(IHIGT(IFIN))-XDONT(IWRK1))
!
!  One takes values <= pivot to ILOWT
!  Again, 2 parts, one where we take care of the remaining
!  high values because we might still need them, and the
!  other when we know that we will have more than enough
!  low values in the end.
!
               JHIG = 0
               do ICRS = 2, IFIN
                  if (XDONT(IHIGT(ICRS)) <= XPIV) then
                     JLOW = JLOW + 1
                     ILOWT (JLOW) = IHIGT (ICRS)
                     if (JLOW >= INTH) exit
                  else
                     JHIG = JHIG + 1
                     IHIGT (JHIG) = IHIGT (ICRS)
                  end if
               end do
!
               do ICRS = ICRS + 1, IFIN
                  if (XDONT(IHIGT(ICRS)) <= XPIV) then
                     JLOW = JLOW + 1
                     ILOWT (JLOW) = IHIGT (ICRS)
                  end if
               end do
            end select
!
!
         case (1)
!
!  Only 1 value is missing in low part
!
            XMIN = XDONT (IHIGT(1))
            IHIG = 1
            do ICRS = 2, JHIG
               if (XDONT(IHIGT(ICRS)) < XMIN) then
                  XMIN = XDONT (IHIGT(ICRS))
                  IHIG = ICRS
               end if
            end do
!
            INDNTH = IHIGT (IHIG)
            return
!
!
         case (0)
!
!  Low part is exactly what we want
!
            exit
!
!
         case (-5:-1)
!
!  Only few values too many in low part
!
            IRNGT (1) = ILOWT (1)
            ILOW = 1 + INTH - JLOW
            do ICRS = 2, INTH
               IWRK = ILOWT (ICRS)
               XWRK = XDONT (IWRK)
               do IDCR = ICRS - 1, MAX (1, ILOW), - 1
                  if (XWRK < XDONT(IRNGT(IDCR))) then
                     IRNGT (IDCR+1) = IRNGT (IDCR)
                  else
                     exit
                  end if
               end do
               IRNGT (IDCR+1) = IWRK
               ILOW = ILOW + 1
            end do
!
            XWRK1 = XDONT (IRNGT(INTH))
            ILOW = 2*INTH - JLOW
            do ICRS = INTH + 1, JLOW
               if (XDONT(ILOWT (ICRS)) < XWRK1) then
                  XWRK = XDONT (ILOWT (ICRS))
                  do IDCR = INTH - 1, MAX (1, ILOW), - 1
                     if (XWRK >= XDONT(IRNGT(IDCR))) exit
                     IRNGT (IDCR+1) = IRNGT (IDCR)
                  end do
                  IRNGT (IDCR+1) = ILOWT (ICRS)
                  XWRK1 = XDONT (IRNGT(INTH))
               end if
               ILOW = ILOW + 1
            end do
!
            INDNTH = IRNGT(INTH)
            return
!
!
         case (:-6)
!
! last case: too many values in low part
!

            IMIL = (JLOW+1) / 2
            IFIN = JLOW
!
!  One chooses a pivot from 1st, last, and middle values
!
            if (XDONT(ILOWT(IMIL)) < XDONT(ILOWT(1))) then
               IWRK = ILOWT (1)
               ILOWT (1) = ILOWT (IMIL)
               ILOWT (IMIL) = IWRK
            end if
            if (XDONT(ILOWT(IMIL)) > XDONT(ILOWT(IFIN))) then
               IWRK = ILOWT (IFIN)
               ILOWT (IFIN) = ILOWT (IMIL)
               ILOWT (IMIL) = IWRK
               if (XDONT(ILOWT(IMIL)) < XDONT(ILOWT(1))) then
                  IWRK = ILOWT (1)
                  ILOWT (1) = ILOWT (IMIL)
                  ILOWT (IMIL) = IWRK
               end if
            end if
            if (IFIN <= 3) exit
!
            XPIV = XDONT (ILOWT(1)) + REAL(INTH)/REAL(JLOW+INTH) * &
                                      (XDONT(ILOWT(IFIN))-XDONT(ILOWT(1)))

!
!  One takes values > XPIV to IHIGT
!
            JHIG = 0
            JLOW = 0
!
            if (XDONT(ILOWT(IFIN)) > XPIV) then
               ICRS = 0
               do
                  ICRS = ICRS + 1
                  if (XDONT(ILOWT(ICRS)) > XPIV) then
                     JHIG = JHIG + 1
                     IHIGT (JHIG) = ILOWT (ICRS)
                     if (ICRS >= IFIN) exit
                  else
                     JLOW = JLOW + 1
                     ILOWT (JLOW) = ILOWT (ICRS)
                     if (JLOW >= INTH) exit
                  end if
               end do
!
               if (ICRS < IFIN) then
                  do
                     ICRS = ICRS + 1
                     if (XDONT(ILOWT(ICRS)) <= XPIV) then
                        JLOW = JLOW + 1
                        ILOWT (JLOW) = ILOWT (ICRS)
                     else
                        if (ICRS >= IFIN) exit
                     end if
                  end do
               end if
            else
               do ICRS = 1, IFIN
                  if (XDONT(ILOWT(ICRS)) > XPIV) then
                     JHIG = JHIG + 1
                     IHIGT (JHIG) = ILOWT (ICRS)
                  else
                     JLOW = JLOW + 1
                     ILOWT (JLOW) = ILOWT (ICRS)
                     if (JLOW >= INTH) exit
                  end if
               end do
!
               do ICRS = ICRS + 1, IFIN
                  if (XDONT(ILOWT(ICRS)) <= XPIV) then
                     JLOW = JLOW + 1
                     ILOWT (JLOW) = ILOWT (ICRS)
                  end if
               end do
            end if
!
         end select
!
      end do
!
!  Now, we only need to find maximum of the 1:INTH set
!

      IWRK1 = ILOWT (1)
      XWRK1 =  XDONT (IWRK1)
      do ICRS = 1+1, INTH
         IWRK = ILOWT (ICRS)
         XWRK = XDONT (IWRK)
         if (XWRK > XWRK1) then
            XWRK1 = XWRK
            IWRK1 = IWRK
         end if
      end do
      INDNTH = IWRK1
      return
!
!
end function D_indnth

pure function R_indnth (XDONT, NORD) result (INDNTH)
!  Return NORDth value of XDONT, i.e fractile of order NORD/SIZE(XDONT).
! __________________________________________________________
!  This routine uses a pivoting strategy such as the one of
!  finding the median based on the quicksort algorithm, but
!  we skew the pivot choice to try to bring it to NORD as
!  fast as possible. It uses 2 temporary arrays, where it
!  stores the indices of the values smaller than the pivot
!  (ILOWT), and the indices of values larger than the pivot
!  that we might still need later on (IHIGT). It iterates
!  until it can bring the number of values in ILOWT to
!  exactly NORD, and then finds the maximum of this set.
!  Michel Olagnon - Aug. 2000
! __________________________________________________________
! _________________________________________________________
      real, dimension (:), intent (in) :: XDONT
      integer :: INDNTH
      integer, intent (in) :: NORD
! __________________________________________________________
      real :: XPIV, XPIV0, XWRK, XWRK1, XMIN, XMAX
!
      integer, dimension (NORD) :: IRNGT
      integer, dimension (SIZE(XDONT)) :: ILOWT, IHIGT
      integer :: NDON, JHIG, JLOW, IHIG, IWRK, IWRK1, IWRK2, IWRK3
      integer :: IMIL, IFIN, ICRS, IDCR, ILOW
      integer :: JLM2, JLM1, JHM2, JHM1, INTH
!
      NDON = SIZE (XDONT)
      INTH = NORD
!
!    First loop is used to fill-in ILOWT, IHIGT at the same time
!
      if (NDON < 2) then
         if (INTH == 1) INDNTH = 1
         return
      end if
!
!  One chooses a pivot, best estimate possible to put fractile near
!  mid-point of the set of low values.
!
      if (XDONT(2) < XDONT(1)) then
         ILOWT (1) = 2
         IHIGT (1) = 1
      else
         ILOWT (1) = 1
         IHIGT (1) = 2
      end if
!
      if (NDON < 3) then
         if (INTH == 1) INDNTH = ILOWT (1)
         if (INTH == 2) INDNTH = IHIGT (1)
         return
      end if
!
      if (XDONT(3) < XDONT(IHIGT(1))) then
         IHIGT (2) = IHIGT (1)
         if (XDONT(3) < XDONT(ILOWT(1))) then
            IHIGT (1) = ILOWT (1)
            ILOWT (1) = 3
         else
            IHIGT (1) = 3
         end if
      else
         IHIGT (2) = 3
      end if
!
      if (NDON < 4) then
         if (INTH == 1) INDNTH = ILOWT (1)
         if (INTH == 2) INDNTH = IHIGT (1)
         if (INTH == 3) INDNTH = IHIGT (2)
         return
      end if
!
      if (XDONT(NDON) < XDONT(IHIGT(1))) then
         IHIGT (3) = IHIGT (2)
         IHIGT (2) = IHIGT (1)
         if (XDONT(NDON) < XDONT(ILOWT(1))) then
            IHIGT (1) = ILOWT (1)
            ILOWT (1) = NDON
         else
            IHIGT (1) = NDON
         end if
      else
         IHIGT (3) = NDON
      end if
!
      if (NDON < 5) then
         if (INTH == 1) INDNTH = ILOWT (1)
         if (INTH == 2) INDNTH = IHIGT (1)
         if (INTH == 3) INDNTH = IHIGT (2)
         if (INTH == 4) INDNTH = IHIGT (3)
         return
      end if
!

      JLOW = 1
      JHIG = 3
      XPIV = XDONT (ILOWT(1)) + REAL(2*INTH)/REAL(NDON+INTH) * &
                                   (XDONT(IHIGT(3))-XDONT(ILOWT(1)))
      if (XPIV >= XDONT(IHIGT(1))) then
         XPIV = XDONT (ILOWT(1)) + REAL(2*INTH)/REAL(NDON+INTH) * &
                                      (XDONT(IHIGT(2))-XDONT(ILOWT(1)))
         if (XPIV >= XDONT(IHIGT(1))) &
             XPIV = XDONT (ILOWT(1)) + REAL (2*INTH) / REAL (NDON+INTH) * &
                                          (XDONT(IHIGT(1))-XDONT(ILOWT(1)))
      end if
      XPIV0 = XPIV
!
!  One puts values > pivot in the end and those <= pivot
!  at the beginning. This is split in 2 cases, so that
!  we can skip the loop test a number of times.
!  As we are also filling in the work arrays at the same time
!  we stop filling in the IHIGT array as soon as we have more
!  than enough values in ILOWT.
!
!
      if (XDONT(NDON) > XPIV) then
         ICRS = 3
         do
            ICRS = ICRS + 1
            if (XDONT(ICRS) > XPIV) then
               if (ICRS >= NDON) exit
               JHIG = JHIG + 1
               IHIGT (JHIG) = ICRS
            else
               JLOW = JLOW + 1
               ILOWT (JLOW) = ICRS
               if (JLOW >= INTH) exit
            end if
         end do
!
!  One restricts further processing because it is no use
!  to store more high values
!
         if (ICRS < NDON-1) then
            do
               ICRS = ICRS + 1
               if (XDONT(ICRS) <= XPIV) then
                  JLOW = JLOW + 1
                  ILOWT (JLOW) = ICRS
               else if (ICRS >= NDON) then
                  exit
               end if
            end do
         end if
!
!
      else
!
!  Same as above, but this is not as easy to optimize, so the
!  DO-loop is kept
!
         do ICRS = 4, NDON - 1
            if (XDONT(ICRS) > XPIV) then
               JHIG = JHIG + 1
               IHIGT (JHIG) = ICRS
            else
               JLOW = JLOW + 1
               ILOWT (JLOW) = ICRS
               if (JLOW >= INTH) exit
            end if
         end do
!
         if (ICRS < NDON-1) then
            do
               ICRS = ICRS + 1
               if (XDONT(ICRS) <= XPIV) then
                  if (ICRS >= NDON) exit
                  JLOW = JLOW + 1
                  ILOWT (JLOW) = ICRS
               end if
            end do
         end if
      end if
!
      JLM2 = 0
      JLM1 = 0
      JHM2 = 0
      JHM1 = 0
      do
         if (JLM2 == JLOW .and. JHM2 == JHIG) then
!
!   We are oscillating. Perturbate by bringing JLOW closer by one
!   to INTH
!
             if (INTH > JLOW) then
                XMIN = XDONT (IHIGT(1))
                IHIG = 1
                do ICRS = 2, JHIG
                   if (XDONT(IHIGT(ICRS)) < XMIN) then
                      XMIN = XDONT (IHIGT(ICRS))
                      IHIG = ICRS
                   end if
                end do
!
                JLOW = JLOW + 1
                ILOWT (JLOW) = IHIGT (IHIG)
                IHIGT (IHIG) = IHIGT (JHIG)
                JHIG = JHIG - 1
             else

                ILOW = ILOWT (1)
                XMAX = XDONT (ILOW)
                do ICRS = 2, JLOW
                   if (XDONT(ILOWT(ICRS)) > XMAX) then
                      IWRK = ILOWT (ICRS)
                      XMAX = XDONT (IWRK)
                      ILOWT (ICRS) = ILOW
                      ILOW = IWRK
                   end if
                end do
                JLOW = JLOW - 1
             end if
         end if
         JLM2 = JLM1
         JLM1 = JLOW
         JHM2 = JHM1
         JHM1 = JHIG
!
!   We try to bring the number of values in the low values set
!   closer to INTH.
!
         select case (INTH-JLOW)
         case (2:)
!
!   Not enough values in low part, at least 2 are missing
!
            INTH = INTH - JLOW
            JLOW = 0
            select case (JHIG)
!!!!!           CASE DEFAULT
!!!!!              write (unit=*,fmt=*) "Assertion failed"
!!!!!              STOP
!
!   We make a special case when we have so few values in
!   the high values set that it is bad performance to choose a pivot
!   and apply the general algorithm.
!
            case (2)
               if (XDONT(IHIGT(1)) <= XDONT(IHIGT(2))) then
                  JLOW = JLOW + 1
                  ILOWT (JLOW) = IHIGT (1)
                  JLOW = JLOW + 1
                  ILOWT (JLOW) = IHIGT (2)
               else
                  JLOW = JLOW + 1
                  ILOWT (JLOW) = IHIGT (2)
                  JLOW = JLOW + 1
                  ILOWT (JLOW) = IHIGT (1)
               end if
               exit
!
            case (3)
!
!
               IWRK1 = IHIGT (1)
               IWRK2 = IHIGT (2)
               IWRK3 = IHIGT (3)
               if (XDONT(IWRK2) < XDONT(IWRK1)) then
                  IHIGT (1) = IWRK2
                  IHIGT (2) = IWRK1
                  IWRK2 = IWRK1
               end if
               if (XDONT(IWRK2) > XDONT(IWRK3)) then
                  IHIGT (3) = IWRK2
                  IHIGT (2) = IWRK3
                  IWRK2 = IWRK3
                  if (XDONT(IWRK2) < XDONT(IHIGT(1))) then
                     IHIGT (2) = IHIGT (1)
                     IHIGT (1) = IWRK2
                  end if
               end if
               JHIG = 0
               do ICRS = JLOW + 1, INTH
                  JHIG = JHIG + 1
                  ILOWT (ICRS) = IHIGT (JHIG)
               end do
               JLOW = INTH
               exit
!
            case (4:)
!
!
               XPIV0 = XPIV
               IFIN = JHIG
!
!  One chooses a pivot from the 2 first values and the last one.
!  This should ensure sufficient renewal between iterations to
!  avoid worst case behavior effects.
!
               IWRK1 = IHIGT (1)
               IWRK2 = IHIGT (2)
               IWRK3 = IHIGT (IFIN)
               if (XDONT(IWRK2) < XDONT(IWRK1)) then
                  IHIGT (1) = IWRK2
                  IHIGT (2) = IWRK1
                  IWRK2 = IWRK1
               end if
               if (XDONT(IWRK2) > XDONT(IWRK3)) then
                  IHIGT (IFIN) = IWRK2
                  IHIGT (2) = IWRK3
                  IWRK2 = IWRK3
                  if (XDONT(IWRK2) < XDONT(IHIGT(1))) then
                     IHIGT (2) = IHIGT (1)
                     IHIGT (1) = IWRK2
                  end if
               end if
!
               IWRK1 = IHIGT (1)
               JLOW = JLOW + 1
               ILOWT (JLOW) = IWRK1
               XPIV = XDONT (IWRK1) + 0.5 * (XDONT(IHIGT(IFIN))-XDONT(IWRK1))
!
!  One takes values <= pivot to ILOWT
!  Again, 2 parts, one where we take care of the remaining
!  high values because we might still need them, and the
!  other when we know that we will have more than enough
!  low values in the end.
!
               JHIG = 0
               do ICRS = 2, IFIN
                  if (XDONT(IHIGT(ICRS)) <= XPIV) then
                     JLOW = JLOW + 1
                     ILOWT (JLOW) = IHIGT (ICRS)
                     if (JLOW >= INTH) exit
                  else
                     JHIG = JHIG + 1
                     IHIGT (JHIG) = IHIGT (ICRS)
                  end if
               end do
!
               do ICRS = ICRS + 1, IFIN
                  if (XDONT(IHIGT(ICRS)) <= XPIV) then
                     JLOW = JLOW + 1
                     ILOWT (JLOW) = IHIGT (ICRS)
                  end if
               end do
            end select
!
!
         case (1)
!
!  Only 1 value is missing in low part
!
            XMIN = XDONT (IHIGT(1))
            IHIG = 1
            do ICRS = 2, JHIG
               if (XDONT(IHIGT(ICRS)) < XMIN) then
                  XMIN = XDONT (IHIGT(ICRS))
                  IHIG = ICRS
               end if
            end do
!
            INDNTH = IHIGT (IHIG)
            return
!
!
         case (0)
!
!  Low part is exactly what we want
!
            exit
!
!
         case (-5:-1)
!
!  Only few values too many in low part
!
            IRNGT (1) = ILOWT (1)
            ILOW = 1 + INTH - JLOW
            do ICRS = 2, INTH
               IWRK = ILOWT (ICRS)
               XWRK = XDONT (IWRK)
               do IDCR = ICRS - 1, MAX (1, ILOW), - 1
                  if (XWRK < XDONT(IRNGT(IDCR))) then
                     IRNGT (IDCR+1) = IRNGT (IDCR)
                  else
                     exit
                  end if
               end do
               IRNGT (IDCR+1) = IWRK
               ILOW = ILOW + 1
            end do
!
            XWRK1 = XDONT (IRNGT(INTH))
            ILOW = 2*INTH - JLOW
            do ICRS = INTH + 1, JLOW
               if (XDONT(ILOWT (ICRS)) < XWRK1) then
                  XWRK = XDONT (ILOWT (ICRS))
                  do IDCR = INTH - 1, MAX (1, ILOW), - 1
                     if (XWRK >= XDONT(IRNGT(IDCR))) exit
                     IRNGT (IDCR+1) = IRNGT (IDCR)
                  end do
                  IRNGT (IDCR+1) = ILOWT (ICRS)
                  XWRK1 = XDONT (IRNGT(INTH))
               end if
               ILOW = ILOW + 1
            end do
!
            INDNTH = IRNGT(INTH)
            return
!
!
         case (:-6)
!
! last case: too many values in low part
!

            IMIL = (JLOW+1) / 2
            IFIN = JLOW
!
!  One chooses a pivot from 1st, last, and middle values
!
            if (XDONT(ILOWT(IMIL)) < XDONT(ILOWT(1))) then
               IWRK = ILOWT (1)
               ILOWT (1) = ILOWT (IMIL)
               ILOWT (IMIL) = IWRK
            end if
            if (XDONT(ILOWT(IMIL)) > XDONT(ILOWT(IFIN))) then
               IWRK = ILOWT (IFIN)
               ILOWT (IFIN) = ILOWT (IMIL)
               ILOWT (IMIL) = IWRK
               if (XDONT(ILOWT(IMIL)) < XDONT(ILOWT(1))) then
                  IWRK = ILOWT (1)
                  ILOWT (1) = ILOWT (IMIL)
                  ILOWT (IMIL) = IWRK
               end if
            end if
            if (IFIN <= 3) exit
!
            XPIV = XDONT (ILOWT(1)) + REAL(INTH)/REAL(JLOW+INTH) * &
                                      (XDONT(ILOWT(IFIN))-XDONT(ILOWT(1)))

!
!  One takes values > XPIV to IHIGT
!
            JHIG = 0
            JLOW = 0
!
            if (XDONT(ILOWT(IFIN)) > XPIV) then
               ICRS = 0
               do
                  ICRS = ICRS + 1
                  if (XDONT(ILOWT(ICRS)) > XPIV) then
                     JHIG = JHIG + 1
                     IHIGT (JHIG) = ILOWT (ICRS)
                     if (ICRS >= IFIN) exit
                  else
                     JLOW = JLOW + 1
                     ILOWT (JLOW) = ILOWT (ICRS)
                     if (JLOW >= INTH) exit
                  end if
               end do
!
               if (ICRS < IFIN) then
                  do
                     ICRS = ICRS + 1
                     if (XDONT(ILOWT(ICRS)) <= XPIV) then
                        JLOW = JLOW + 1
                        ILOWT (JLOW) = ILOWT (ICRS)
                     else
                        if (ICRS >= IFIN) exit
                     end if
                  end do
               end if
            else
               do ICRS = 1, IFIN
                  if (XDONT(ILOWT(ICRS)) > XPIV) then
                     JHIG = JHIG + 1
                     IHIGT (JHIG) = ILOWT (ICRS)
                  else
                     JLOW = JLOW + 1
                     ILOWT (JLOW) = ILOWT (ICRS)
                     if (JLOW >= INTH) exit
                  end if
               end do
!
               do ICRS = ICRS + 1, IFIN
                  if (XDONT(ILOWT(ICRS)) <= XPIV) then
                     JLOW = JLOW + 1
                     ILOWT (JLOW) = ILOWT (ICRS)
                  end if
               end do
            end if
!
         end select
!
      end do
!
!  Now, we only need to find maximum of the 1:INTH set
!

      IWRK1 = ILOWT (1)
      XWRK1 =  XDONT (IWRK1)
      do ICRS = 1+1, INTH
         IWRK = ILOWT (ICRS)
         XWRK = XDONT (IWRK)
         if (XWRK > XWRK1) then
            XWRK1 = XWRK
            IWRK1 = IWRK
         end if
      end do
      INDNTH = IWRK1
      return
!
!
end function R_indnth
pure function I_indnth (XDONT, NORD) result (INDNTH)
!  Return NORDth value of XDONT, i.e fractile of order NORD/SIZE(XDONT).
! __________________________________________________________
!  This routine uses a pivoting strategy such as the one of
!  finding the median based on the quicksort algorithm, but
!  we skew the pivot choice to try to bring it to NORD as
!  fast as possible. It uses 2 temporary arrays, where it
!  stores the indices of the values smaller than the pivot
!  (ILOWT), and the indices of values larger than the pivot
!  that we might still need later on (IHIGT). It iterates
!  until it can bring the number of values in ILOWT to
!  exactly NORD, and then finds the maximum of this set.
!  Michel Olagnon - Aug. 2000
! __________________________________________________________
! __________________________________________________________
      integer, dimension (:), intent (in)  :: XDONT
      integer :: INDNTH
      integer, intent (in) :: NORD
! __________________________________________________________
      integer :: XPIV, XPIV0, XWRK, XWRK1, XMIN, XMAX
!
      integer, dimension (NORD) :: IRNGT
      integer, dimension (SIZE(XDONT)) :: ILOWT, IHIGT
      integer :: NDON, JHIG, JLOW, IHIG, IWRK, IWRK1, IWRK2, IWRK3
      integer :: IMIL, IFIN, ICRS, IDCR, ILOW
      integer :: JLM2, JLM1, JHM2, JHM1, INTH
!
      NDON = SIZE (XDONT)
      INTH = NORD
!
!    First loop is used to fill-in ILOWT, IHIGT at the same time
!
      if (NDON < 2) then
         if (INTH == 1) INDNTH = 1
         return
      end if
!
!  One chooses a pivot, best estimate possible to put fractile near
!  mid-point of the set of low values.
!
      if (XDONT(2) < XDONT(1)) then
         ILOWT (1) = 2
         IHIGT (1) = 1
      else
         ILOWT (1) = 1
         IHIGT (1) = 2
      end if
!
      if (NDON < 3) then
         if (INTH == 1) INDNTH = ILOWT (1)
         if (INTH == 2) INDNTH = IHIGT (1)
         return
      end if
!
      if (XDONT(3) < XDONT(IHIGT(1))) then
         IHIGT (2) = IHIGT (1)
         if (XDONT(3) < XDONT(ILOWT(1))) then
            IHIGT (1) = ILOWT (1)
            ILOWT (1) = 3
         else
            IHIGT (1) = 3
         end if
      else
         IHIGT (2) = 3
      end if
!
      if (NDON < 4) then
         if (INTH == 1) INDNTH = ILOWT (1)
         if (INTH == 2) INDNTH = IHIGT (1)
         if (INTH == 3) INDNTH = IHIGT (2)
         return
      end if
!
      if (XDONT(NDON) < XDONT(IHIGT(1))) then
         IHIGT (3) = IHIGT (2)
         IHIGT (2) = IHIGT (1)
         if (XDONT(NDON) < XDONT(ILOWT(1))) then
            IHIGT (1) = ILOWT (1)
            ILOWT (1) = NDON
         else
            IHIGT (1) = NDON
         end if
      else
         IHIGT (3) = NDON
      end if
!
      if (NDON < 5) then
         if (INTH == 1) INDNTH = ILOWT (1)
         if (INTH == 2) INDNTH = IHIGT (1)
         if (INTH == 3) INDNTH = IHIGT (2)
         if (INTH == 4) INDNTH = IHIGT (3)
         return
      end if
!

      JLOW = 1
      JHIG = 3
      XPIV = XDONT (ILOWT(1)) + REAL(2*INTH)/REAL(NDON+INTH) * &
                                   (XDONT(IHIGT(3))-XDONT(ILOWT(1)))
      if (XPIV >= XDONT(IHIGT(1))) then
         XPIV = XDONT (ILOWT(1)) + REAL(2*INTH)/REAL(NDON+INTH) * &
                                      (XDONT(IHIGT(2))-XDONT(ILOWT(1)))
         if (XPIV >= XDONT(IHIGT(1))) &
             XPIV = XDONT (ILOWT(1)) + REAL (2*INTH) / REAL (NDON+INTH) * &
                                          (XDONT(IHIGT(1))-XDONT(ILOWT(1)))
      end if
      XPIV0 = XPIV
!
!  One puts values > pivot in the end and those <= pivot
!  at the beginning. This is split in 2 cases, so that
!  we can skip the loop test a number of times.
!  As we are also filling in the work arrays at the same time
!  we stop filling in the IHIGT array as soon as we have more
!  than enough values in ILOWT.
!
!
      if (XDONT(NDON) > XPIV) then
         ICRS = 3
         do
            ICRS = ICRS + 1
            if (XDONT(ICRS) > XPIV) then
               if (ICRS >= NDON) exit
               JHIG = JHIG + 1
               IHIGT (JHIG) = ICRS
            else
               JLOW = JLOW + 1
               ILOWT (JLOW) = ICRS
               if (JLOW >= INTH) exit
            end if
         end do
!
!  One restricts further processing because it is no use
!  to store more high values
!
         if (ICRS < NDON-1) then
            do
               ICRS = ICRS + 1
               if (XDONT(ICRS) <= XPIV) then
                  JLOW = JLOW + 1
                  ILOWT (JLOW) = ICRS
               else if (ICRS >= NDON) then
                  exit
               end if
            end do
         end if
!
!
      else
!
!  Same as above, but this is not as easy to optimize, so the
!  DO-loop is kept
!
         do ICRS = 4, NDON - 1
            if (XDONT(ICRS) > XPIV) then
               JHIG = JHIG + 1
               IHIGT (JHIG) = ICRS
            else
               JLOW = JLOW + 1
               ILOWT (JLOW) = ICRS
               if (JLOW >= INTH) exit
            end if
         end do
!
         if (ICRS < NDON-1) then
            do
               ICRS = ICRS + 1
               if (XDONT(ICRS) <= XPIV) then
                  if (ICRS >= NDON) exit
                  JLOW = JLOW + 1
                  ILOWT (JLOW) = ICRS
               end if
            end do
         end if
      end if
!
      JLM2 = 0
      JLM1 = 0
      JHM2 = 0
      JHM1 = 0
      do
         if (JLM2 == JLOW .and. JHM2 == JHIG) then
!
!   We are oscillating. Perturbate by bringing JLOW closer by one
!   to INTH
!
             if (INTH > JLOW) then
                XMIN = XDONT (IHIGT(1))
                IHIG = 1
                do ICRS = 2, JHIG
                   if (XDONT(IHIGT(ICRS)) < XMIN) then
                      XMIN = XDONT (IHIGT(ICRS))
                      IHIG = ICRS
                   end if
                end do
!
                JLOW = JLOW + 1
                ILOWT (JLOW) = IHIGT (IHIG)
                IHIGT (IHIG) = IHIGT (JHIG)
                JHIG = JHIG - 1
             else

                ILOW = ILOWT (1)
                XMAX = XDONT (ILOW)
                do ICRS = 2, JLOW
                   if (XDONT(ILOWT(ICRS)) > XMAX) then
                      IWRK = ILOWT (ICRS)
                      XMAX = XDONT (IWRK)
                      ILOWT (ICRS) = ILOW
                      ILOW = IWRK
                   end if
                end do
                JLOW = JLOW - 1
             end if
         end if
         JLM2 = JLM1
         JLM1 = JLOW
         JHM2 = JHM1
         JHM1 = JHIG
!
!   We try to bring the number of values in the low values set
!   closer to INTH.
!
         select case (INTH-JLOW)
         case (2:)
!
!   Not enough values in low part, at least 2 are missing
!
            INTH = INTH - JLOW
            JLOW = 0
            select case (JHIG)
!!!!!           CASE DEFAULT
!!!!!              write (unit=*,fmt=*) "Assertion failed"
!!!!!              STOP
!
!   We make a special case when we have so few values in
!   the high values set that it is bad performance to choose a pivot
!   and apply the general algorithm.
!
            case (2)
               if (XDONT(IHIGT(1)) <= XDONT(IHIGT(2))) then
                  JLOW = JLOW + 1
                  ILOWT (JLOW) = IHIGT (1)
                  JLOW = JLOW + 1
                  ILOWT (JLOW) = IHIGT (2)
               else
                  JLOW = JLOW + 1
                  ILOWT (JLOW) = IHIGT (2)
                  JLOW = JLOW + 1
                  ILOWT (JLOW) = IHIGT (1)
               end if
               exit
!
            case (3)
!
!
               IWRK1 = IHIGT (1)
               IWRK2 = IHIGT (2)
               IWRK3 = IHIGT (3)
               if (XDONT(IWRK2) < XDONT(IWRK1)) then
                  IHIGT (1) = IWRK2
                  IHIGT (2) = IWRK1
                  IWRK2 = IWRK1
               end if
               if (XDONT(IWRK2) > XDONT(IWRK3)) then
                  IHIGT (3) = IWRK2
                  IHIGT (2) = IWRK3
                  IWRK2 = IWRK3
                  if (XDONT(IWRK2) < XDONT(IHIGT(1))) then
                     IHIGT (2) = IHIGT (1)
                     IHIGT (1) = IWRK2
                  end if
               end if
               JHIG = 0
               do ICRS = JLOW + 1, INTH
                  JHIG = JHIG + 1
                  ILOWT (ICRS) = IHIGT (JHIG)
               end do
               JLOW = INTH
               exit
!
            case (4:)
!
!
               XPIV0 = XPIV
               IFIN = JHIG
!
!  One chooses a pivot from the 2 first values and the last one.
!  This should ensure sufficient renewal between iterations to
!  avoid worst case behavior effects.
!
               IWRK1 = IHIGT (1)
               IWRK2 = IHIGT (2)
               IWRK3 = IHIGT (IFIN)
               if (XDONT(IWRK2) < XDONT(IWRK1)) then
                  IHIGT (1) = IWRK2
                  IHIGT (2) = IWRK1
                  IWRK2 = IWRK1
               end if
               if (XDONT(IWRK2) > XDONT(IWRK3)) then
                  IHIGT (IFIN) = IWRK2
                  IHIGT (2) = IWRK3
                  IWRK2 = IWRK3
                  if (XDONT(IWRK2) < XDONT(IHIGT(1))) then
                     IHIGT (2) = IHIGT (1)
                     IHIGT (1) = IWRK2
                  end if
               end if
!
               IWRK1 = IHIGT (1)
               JLOW = JLOW + 1
               ILOWT (JLOW) = IWRK1
               XPIV = XDONT (IWRK1) + 0.5 * (XDONT(IHIGT(IFIN))-XDONT(IWRK1))
!
!  One takes values <= pivot to ILOWT
!  Again, 2 parts, one where we take care of the remaining
!  high values because we might still need them, and the
!  other when we know that we will have more than enough
!  low values in the end.
!
               JHIG = 0
               do ICRS = 2, IFIN
                  if (XDONT(IHIGT(ICRS)) <= XPIV) then
                     JLOW = JLOW + 1
                     ILOWT (JLOW) = IHIGT (ICRS)
                     if (JLOW >= INTH) exit
                  else
                     JHIG = JHIG + 1
                     IHIGT (JHIG) = IHIGT (ICRS)
                  end if
               end do
!
               do ICRS = ICRS + 1, IFIN
                  if (XDONT(IHIGT(ICRS)) <= XPIV) then
                     JLOW = JLOW + 1
                     ILOWT (JLOW) = IHIGT (ICRS)
                  end if
               end do
            end select
!
!
         case (1)
!
!  Only 1 value is missing in low part
!
            XMIN = XDONT (IHIGT(1))
            IHIG = 1
            do ICRS = 2, JHIG
               if (XDONT(IHIGT(ICRS)) < XMIN) then
                  XMIN = XDONT (IHIGT(ICRS))
                  IHIG = ICRS
               end if
            end do
!
            INDNTH = IHIGT (IHIG)
            return
!
!
         case (0)
!
!  Low part is exactly what we want
!
            exit
!
!
         case (-5:-1)
!
!  Only few values too many in low part
!
            IRNGT (1) = ILOWT (1)
            ILOW = 1 + INTH - JLOW
            do ICRS = 2, INTH
               IWRK = ILOWT (ICRS)
               XWRK = XDONT (IWRK)
               do IDCR = ICRS - 1, MAX (1, ILOW), - 1
                  if (XWRK < XDONT(IRNGT(IDCR))) then
                     IRNGT (IDCR+1) = IRNGT (IDCR)
                  else
                     exit
                  end if
               end do
               IRNGT (IDCR+1) = IWRK
               ILOW = ILOW + 1
            end do
!
            XWRK1 = XDONT (IRNGT(INTH))
            ILOW = 2*INTH - JLOW
            do ICRS = INTH + 1, JLOW
               if (XDONT(ILOWT (ICRS)) < XWRK1) then
                  XWRK = XDONT (ILOWT (ICRS))
                  do IDCR = INTH - 1, MAX (1, ILOW), - 1
                     if (XWRK >= XDONT(IRNGT(IDCR))) exit
                     IRNGT (IDCR+1) = IRNGT (IDCR)
                  end do
                  IRNGT (IDCR+1) = ILOWT (ICRS)
                  XWRK1 = XDONT (IRNGT(INTH))
               end if
               ILOW = ILOW + 1
            end do
!
            INDNTH = IRNGT(INTH)
            return
!
!
         case (:-6)
!
! last case: too many values in low part
!

            IMIL = (JLOW+1) / 2
            IFIN = JLOW
!
!  One chooses a pivot from 1st, last, and middle values
!
            if (XDONT(ILOWT(IMIL)) < XDONT(ILOWT(1))) then
               IWRK = ILOWT (1)
               ILOWT (1) = ILOWT (IMIL)
               ILOWT (IMIL) = IWRK
            end if
            if (XDONT(ILOWT(IMIL)) > XDONT(ILOWT(IFIN))) then
               IWRK = ILOWT (IFIN)
               ILOWT (IFIN) = ILOWT (IMIL)
               ILOWT (IMIL) = IWRK
               if (XDONT(ILOWT(IMIL)) < XDONT(ILOWT(1))) then
                  IWRK = ILOWT (1)
                  ILOWT (1) = ILOWT (IMIL)
                  ILOWT (IMIL) = IWRK
               end if
            end if
            if (IFIN <= 3) exit
!
            XPIV = XDONT (ILOWT(1)) + REAL(INTH)/REAL(JLOW+INTH) * &
                                      (XDONT(ILOWT(IFIN))-XDONT(ILOWT(1)))

!
!  One takes values > XPIV to IHIGT
!
            JHIG = 0
            JLOW = 0
!
            if (XDONT(ILOWT(IFIN)) > XPIV) then
               ICRS = 0
               do
                  ICRS = ICRS + 1
                  if (XDONT(ILOWT(ICRS)) > XPIV) then
                     JHIG = JHIG + 1
                     IHIGT (JHIG) = ILOWT (ICRS)
                     if (ICRS >= IFIN) exit
                  else
                     JLOW = JLOW + 1
                     ILOWT (JLOW) = ILOWT (ICRS)
                     if (JLOW >= INTH) exit
                  end if
               end do
!
               if (ICRS < IFIN) then
                  do
                     ICRS = ICRS + 1
                     if (XDONT(ILOWT(ICRS)) <= XPIV) then
                        JLOW = JLOW + 1
                        ILOWT (JLOW) = ILOWT (ICRS)
                     else
                        if (ICRS >= IFIN) exit
                     end if
                  end do
               end if
            else
               do ICRS = 1, IFIN
                  if (XDONT(ILOWT(ICRS)) > XPIV) then
                     JHIG = JHIG + 1
                     IHIGT (JHIG) = ILOWT (ICRS)
                  else
                     JLOW = JLOW + 1
                     ILOWT (JLOW) = ILOWT (ICRS)
                     if (JLOW >= INTH) exit
                  end if
               end do
!
               do ICRS = ICRS + 1, IFIN
                  if (XDONT(ILOWT(ICRS)) <= XPIV) then
                     JLOW = JLOW + 1
                     ILOWT (JLOW) = ILOWT (ICRS)
                  end if
               end do
            end if
!
         end select
!
      end do
!
!  Now, we only need to find maximum of the 1:INTH set
!

      IWRK1 = ILOWT (1)
      XWRK1 =  XDONT (IWRK1)
      do ICRS = 1+1, INTH
         IWRK = ILOWT (ICRS)
         XWRK = XDONT (IWRK)
         if (XWRK > XWRK1) then
            XWRK1 = XWRK
            IWRK1 = IWRK
         end if
      end do
      INDNTH = IWRK1
      return
!
!
end function I_indnth

end module sorting_partial_mod
