// ENDHEADER

//
// Basic usage is printf-like,
// %f     - floataing point, default precision
//
// % [ flag ] [ width ] [ .precision ] [ conversion ]
//
// Allowed conversion characters are d,e,f
// Terminal colours are also allowed, in the form of an escape code
// \[spec]{text}
// where spec is a comma-separated list (may be empty) of either
// integer codes or string descriptive names (not case sensitive)
// as listed in terminal::color.
// A color spec can start with either '\[' (escape code), or with '\\' '[' sequence

#if !defined(MPTOOLKIT_COMMON_OUTPUT_FORMATTER_H)
#define MPTOOLKIT_COMMON_OUTPUT_FORMATTER_H

enum class AllowColor { Never, Auto, Always };

class OutputFormatter
{
   public:


      void ShowColors(AllowColor c);

      template <typename... Args>
      std::string Format(std::string FormatSpec, Args);
};

struct FormatSpec
{
   char Flag;  // or \0 if no flag
   int Width;  // 0 if no width specified
   int Precision; // -1 if no precision specified
   char Conversion;  // 'd', 'e', or 'f' for standard printf, or a generic character

   FormatSpec();

   // precondition: *beg == '%' && *(beg+1) != '%'
   // Parses a FormatSpec at the current location in the string, and
   // returns the FormatSpec object, advances the beg iterator to one-past-the-end
   // of the FormatSpec
   static FormatSpec parse(char const*& beg, char const* end);

   std::string to_string() const;

   static is_flag(char c);

};

FormatSpec 
FormatSpec::parse(char const*& beg, char const* end)
{
   FormatSpec Result;
   ++beg;
   CHECK(beg < end);
   if (is_flag(*beg))
   {
      Result.Fiag = *beg++;
   }
   CHECK(beg < end);
   if (*beg >= '0' && *beg <= '9')
   {
      Result.Width = std::strtoi(beg, &beg, 10);
   }
   CHECK(beg < end);
   if (*beg == '.')
   {
      ++beg;
      Result.Precision = std::strtoi(beg, &beg, 10);
   }
   CHECK(beg < end);
   Result.Conversion = *beg++;
   return Result;
}

std::string
FormatSpec::to_string() const
{
   std::string Out = "%";
   if (Flag != '\0')
      Out += Flag;
   if (Width != 0)
      Out += std::to_string(Width);
   if (Precision != -1)
      Out += std::to_string(Precision);
   Out += Conversion;
   return Out;
}
   
template <typename... Args>
std::string
OutputFormatter::Format(std::string FormatSpec, Args)
{

}

