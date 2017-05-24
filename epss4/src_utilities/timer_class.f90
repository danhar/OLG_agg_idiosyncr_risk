!
!////////////////////////////////////////////////////////////////////////
!      timer_class.f90
!      Created: September 2, 2013 1:58 PM By: NocturnalAviationSoftware  
!
!      See description and example on:
!      http://www.nocturnalaviationsoftware.com/free-fortran-code.html
!
!      Defines a class for timing Fortran program execution
!      Attention: uses cpu_time (which differs depending on how many threads there are) , not system_clock (which measures wall clock time)
!      Therefore good for timing algorithms /code sections, in particular PURE ones.
!
!      Usage
!         * Initialization *
!               CALL timer % init()
!         * Starting the timer *
!               CALL timer % start()
!         * Stopping the timer *
!               CALL timer % stop()
!         * Reading the time *
!               time = timer % elapsedTime(units)
!           units (optional) = TC_SECONDS or TC_MINUTES or TC_HOURS
!
!       Example
!           PROGRAM timeIt 
!               use timer_class ,only: t_timer, tc_minutes
!               type(t_timer) :: innerTimer, outerTimer
!               call outerTimer%start()
!               ! Do some work
!               call innerTimer%start()
!               ! Do some more work
!               call innerTimer%stop()
!               ! Do even more work
!               call outerTimer%stop()
!               print *,"Outer work=",outerTimer%elapsedTime(tc_minutes)," min"
!               print *, "Inner work = ", innerTimer % elapsedTime()," sec"
!           END PROGRAM timeIt
!
!////////////////////////////////////////////////////////////////////////
!
      MODULE timer_class
      IMPLICIT NONE
      PRIVATE
!
!     -----------------
!     Private constants
!     -----------------
!
      INTEGER, PARAMETER, PRIVATE :: d = 15
!
!     ----------------
!     Public constants
!     ----------------
!
      INTEGER, PARAMETER, PUBLIC  :: TP = SELECTED_REAL_KIND(d)
      INTEGER, PARAMETER, PUBLIC  :: TC_SECONDS = 0, TC_MINUTES = 1, TC_HOURS = 2
!
!     ---------------------
!     Class type definition
!     ---------------------
!
      TYPE, PUBLIC :: t_timer
         LOGICAL      , PRIVATE :: started = .FALSE., stopped = .FALSE.
         REAL(KIND=TP), PRIVATE :: startTime = 0.0_TP
         REAL(KIND=TP), PRIVATE :: finishTime = 0.0_TP
!
!        ========
         CONTAINS
!        ========
!
         PROCEDURE, PASS :: start => startTimer
         PROCEDURE, PASS :: stop  => stopTimer
         PROCEDURE, PASS :: elapsedTime
         
      END TYPE t_timer
!
!     ========
      CONTAINS
!     ========
! 
!
!////////////////////////////////////////////////////////////////////////  
!  
      SUBROUTINE startTimer(self)  
         IMPLICIT NONE
         CLASS(t_timer) :: self
         self % started = .TRUE.
         CALL CPU_TIME(self % startTime)         
      END SUBROUTINE startTimer
!
!////////////////////////////////////////////////////////////////////////  
!  
      SUBROUTINE stopTimer(self)  
         IMPLICIT NONE
         CLASS(t_timer) :: self
         CALL CPU_TIME(self % finishTime)
         self % stopped = .TRUE.
      END SUBROUTINE stopTimer
!
!//////////////////////////////////////////////////////////////////////// 
! 
       FUNCTION elapsedTime(self,units)  
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         CLASS(t_timer)    :: self
         INTEGER, OPTIONAL :: units
         REAL(KIND=TP)  :: elapsedTime
!
!        ------------------------------------------
!        Return zero if the timer was never started
!        ------------------------------------------
!
         IF ( .NOT.self % started )     THEN
            elapsedTime = 0.0_TP
            RETURN
         END IF 
!
!        ----------------------------------------------
!        If the timer was not stopped, then return the 
!        current time elapsed
!        ----------------------------------------------
!
         IF ( .NOT.self % stopped )     THEN
            CALL self % stop() 
         END IF 

         elapsedTime =  self % finishTime - self % startTime
!
!        -------------------------------------
!        Convert to requested units if present
!        -------------------------------------
!
         IF ( PRESENT(units) )     THEN
         
            SELECT CASE ( units )
               CASE( TC_MINUTES ) 
                  elapsedTime = elapsedTime/60.0_TP
               CASE( TC_HOURS )
                  elapsedTime = elapsedTime/3600.0_TP
               CASE DEFAULT 
               
            END SELECT 
         END IF 
      
      END FUNCTION elapsedTime
      
      END MODULE timer_class
