#ifndef UTILS_TEMPLATES_H
#define UTILS_TEMPLATES_H

namespace mp {

/// Source - https://stackoverflow.com/a/43526890
/// Posted by R Sahu, modified by community. See post 'Timeline' for change history
/// Retrieved 2025-11-14, License - CC BY-SA 3.0
template <typename> struct FirstArgument;

/// Specialize for function
template <typename R, typename A, typename... Args>
struct FirstArgument<R(A, Args...)>
{
  using type = A;
};

/// Argument type access
template <typename T>
using first_argument_t = typename FirstArgument<T>::type;

/*
void foo(int a){ }

void bar(int a, double b){ }

int main()
{
  long value = 1L;
  foo(static_cast<first_agument_t<decltype(foo)>>(value) );
  bar(static_cast<first_agument_t<decltype(bar)>>(value), 0);
}
*/

}  // namespace mp

#endif // UTILS_TEMPLATES_H
