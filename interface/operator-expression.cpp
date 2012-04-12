

std::pair<OperatorList, MPOperator>
ParseLatticeAndOperator(std::string const& str)
{
   
}

typedef std::complex<double> complex;

typedef boost::variant<complex, MPOperator> ElementType;

struct calculator : grammar<calculator>
{
   
   LatticeFile = 

   real_number

   imag_number

   mp_operator

   factor

   term

   expression =
      term
   >> *( ('+' >> term)[make_op(plus<element_type>(), self.eval)]
         | ('-' >> term)[make_op(minus<element_type>(), self.eval)]
         )
      ;


      
