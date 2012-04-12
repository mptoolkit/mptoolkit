// -*- C++ -*- $Id$

template <typename T>
Lattice::Lattice(SiteBlock const& s, T const& Coord)
   : Data_(s), 
     Coordinates_(1, boost::lexical_cast<std::string>(Coord)), 
     CoordinatesFixed_(false)
{
}

template <typename Visitor>
typename Visitor::result_type
Lattice::apply_visitor(Visitor const& v) const
{
   return Data_.apply_visitor(v);
}

template <typename Visitor>
typename Visitor::result_type
Lattice::apply_visitor(Visitor const& v)
{
   return Data_.apply_visitor(v);
}

template <typename T1>
int Lattice::site_at_coordinate(T1 const& x1) const
{
   return this->site_at_coordinate(boost::lexical_cast<std::string>(x1));
}

template <typename T1, typename T2>
int Lattice::site_at_coordinate(T1 const& x1, T2 const& x2) const
{
   return this->site_at_coordinate(boost::lexical_cast<std::string>(x1)+','
                                   +boost::lexical_cast<std::string>(x2));
}

template <typename T1, typename T2, typename T3>
int Lattice::site_at_coordinate(T1 const& x1, T2 const& x2, T3 const& x3) const
{
   return this->site_at_coordinate(boost::lexical_cast<std::string>(x1)+','
                                   +boost::lexical_cast<std::string>(x2)+','
                                   +boost::lexical_cast<std::string>(x3));
}

template <typename FwdIter>
void Lattice::fix_coordinates_from_sequence(FwdIter start, FwdIter finish)
{
   Coordinates_ = std::vector<std::string>(start, finish);
   CHECK_EQUAL(Data_.size(), int(Coordinates_.size()))
      ("The list of coordinates must be the same size as the lattice");
   this->Fixate();
}

template <typename T>
void Lattice::fix_coordinates(T const& c1)
{
   CHECK_EQUAL(Data_.size(), 1);
   Coordinates_.clear();
   Coordinates_.push_back(boost::lexical_cast<std::string>(c1));
   this->Fixate();
}

template <typename T>
void Lattice::fix_coordinates(T const& c1, T const& c2)
{
   CHECK_EQUAL(Data_.size(), 2);
   Coordinates_.clear();
   Coordinates_.push_back(boost::lexical_cast<std::string>(c1));
   Coordinates_.push_back(boost::lexical_cast<std::string>(c2));
   this->Fixate();
}

template <typename T>
void Lattice::fix_coordinates(T const& c1, T const& c2, T const& c3)
{
   CHECK_EQUAL(Data_.size(), 3);
   Coordinates_.clear();
   Coordinates_.push_back(boost::lexical_cast<std::string>(c1));
   Coordinates_.push_back(boost::lexical_cast<std::string>(c2));
   Coordinates_.push_back(boost::lexical_cast<std::string>(c3));
   this->Fixate();
}

template <typename T>
void Lattice::fix_coordinates(T const& c1, T const& c2, T const& c3, T const& c4)
{
   CHECK_EQUAL(Data_.size(), 4);
   Coordinates_.clear();
   Coordinates_.push_back(boost::lexical_cast<std::string>(c1));
   Coordinates_.push_back(boost::lexical_cast<std::string>(c2));
   Coordinates_.push_back(boost::lexical_cast<std::string>(c3));
   Coordinates_.push_back(boost::lexical_cast<std::string>(c4));
   this->Fixate();
}
