#
#  Copyright (C) 2025 RDKit contributors
#
#  All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#     * Neither the name of the RDKit contributors nor the names of
#       contributors may be used to endorse or promote products derived
#       from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#

import unittest
import os

from rdkit import Chem, RDConfig
from rdkit.Chem import rdChemReactions


class TestRdfParser(unittest.TestCase):

  def setUp(self):
    # Sample RDF block from the C++ tests
    self.sampleRdf = """$RDFILE 1
$DATM 12/9/2025, 1:46:53 PM
$RFMT
$RXN

      RDKit

  2  1
$MOL

     RDKit          2D

 10 10  0  0  0  0  0  0  0  0999 V2000
   -0.2579   -0.8086    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6063   -0.3056    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6027    0.6944    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4741   -0.8024    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.3383   -0.2994    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.2063   -0.7964    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0703   -0.2932    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0669    0.7068    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.1989    1.2038    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.3349    0.7006    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  2  4  1  0
  4  5  1  0
  5  6  2  0
  6  7  1  0
  7  8  2  0
  8  9  1  0
  9 10  2  0
 10  5  1  0
M  END
$MOL

     RDKit          2D

  5  4  0  0  0  0  0  0  0  0999 V2000
    6.6267    1.2014    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.6287    0.2014    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.7637   -0.3004    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.8967    0.1980    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    5.7657   -1.3004    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  3  4  2  0
  3  5  1  0
M  END
$MOL

     RDKit          2D

 15 15  0  0  0  0  0  0  0  0999 V2000
    7.2537   -1.4305    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.1201   -0.9313    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.1209    0.0687    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.9857   -1.4319    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.8519   -0.9327    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   10.7177   -1.4333    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   11.5841   -0.9341    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   11.5849    0.0661    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   12.4513    0.5653    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   13.3169    0.0645    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   12.4521    1.5653    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   11.5865    2.0659    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   13.3185    2.0645    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   10.7193    0.5665    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.8527    0.0673    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  2  4  1  0
  4  5  1  0
  5  6  2  0
  6  7  1  0
  7  8  2  0
  8  9  1  0
  9 10  1  0
  9 11  1  0
 11 12  2  0
 11 13  1  0
  8 14  1  0
 14 15  2  0
 15  5  1  0
M  END
"""

  def test_RdfBlockIsReaction(self):
    """Test that RdfBlockIsReaction correctly identifies RDF blocks"""
    self.assertTrue(rdChemReactions.RdfBlockIsReaction(self.sampleRdf))
    self.assertFalse(rdChemReactions.RdfBlockIsReaction("not an RDF block"))

  def test_ReactionsFromRdfBlock(self):
    """Test parsing reactions from RDF block"""
    rxns = rdChemReactions.ReactionsFromRdfBlock(self.sampleRdf)
    self.assertIsNotNone(rxns)
    self.assertGreater(len(rxns), 0)
    
    # Check that we got a reaction
    rxn = rxns[0]
    self.assertIsNotNone(rxn)
    
    # Check that the reaction has reactants and products
    self.assertGreater(rxn.GetNumReactantTemplates(), 0)
    self.assertGreater(rxn.GetNumProductTemplates(), 0)


if __name__ == '__main__':
  unittest.main()
