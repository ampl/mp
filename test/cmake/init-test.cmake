function(expect_eq expected actual)
  if (NOT "${expected}" STREQUAL "${actual}")
    message(SEND_ERROR "Test failed.\nExpected: '${expected}'\nActual: '${actual}'")
  endif ()
endfunction()

include(init)

join(result)
expect_eq("" "${result}")

join(result ab cd ef)
expect_eq("ab cd ef" "${result}")

set_cache(test a STRING doc string)
expect_eq(a "${test}")

set_cache(test b STRING doc)
expect_eq(a "${test}")

set_cache(test b STRING doc FORCE)
expect_eq(b "${test}")
