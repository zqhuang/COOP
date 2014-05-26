module coop_basictype
  implicit none
!!$private
!!$public:: coop_link, coop_list
!!$
!!$  type coop_link
!!$     private
!!$     class(*),pointer:: value => null()
!!$     type(coop_link),pointer:: next ==> null()
!!$   contains
!!$     procedure::getValue => coop_link_getValue
!!$     procedure::nextLink => coop_link_nextLink
!!$     procedure::setNestLink => coop_link_setNextLink
!!$  end type coop_link
!!$
!!$  interface coop_link
!!$     procedure coop_link_constructor
!!$  end interface coop_link
!!$
!!$
!!$contains
!!$
!!$  function coop_link_nextLink(this) result(nextLink)
!!$    class(coop_link):: this
!!$    class(coop_link),pointer::nextLink
!!$    nextLink => this%next
!!$  end function coop_link_nextLink
!!$  
!!$  subroutine coop_link_setNextLink(this, next)
!!$    class(coop_link) :: this
!!$    class(coop_link), pointer:: next
!!$    this%next => next
!!$  end subroutine coop_link_setNextLink
!!$
!!$  function coop_link_getValue(this) result(getValue)
!!$    class(coop_link):: this
!!$    class(*),pointer:: getValue
!!$    getValue => this%value
!!$  end function coop_link_getValue
!!$
!!$  function coop_link_constructor(value, next) result(constructor)
!!$    class(link),pointer:: constructor
!!$    class(*) :: value
!!$    class(link),pointer::next
!!$    allocate(constructor)
!!$    constructor%next => next
!!$    allocate(constructor%value, source = value)
!!$  end function coop_link_constructor
!!$
!!$
end module coop_basictype
