Important information about the code
==

This folder contains the source of Afivo, split over several modules. The
m_aX_... files contain code for both the 2D and 3D version. The m_a2_... and
m_a3_... files are generated automatically, in the following way:

1. Replace $D by 2 or 3 (dimension of code)
2. Preprocess file with cpp
3. cat -s (merge multiple blank lines)

The m_a2_... and m_a3_... files are automaticcaly generated from the m_aX_...
files. This means that you should **never** edit the m_a2_... or m_a3_... files,
but only the m_aX_... files!
