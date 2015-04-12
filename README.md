# Tulip Gomory-Hu tree plugin

## Description

This plugin allow to compute the [Gomory-Hu tree](http://en.wikipedia.org/wiki/Gomory%E2%80%93Hu_tree) associated to a weighted undirected graph.

This implementation follows the notation available in the [Gomory-Hu tree](http://en.wikipedia.org/wiki/Gomory%E2%80%93Hu_tree) Wikipedia page, and relies on the `boost::boykov_kolmogorov_max_flow` algorithm from the [Boost Graph Library](http://www.boost.org/doc/libs/1_57_0/libs/graph/doc/boykov_kolmogorov_max_flow.html) to compute minimum cuts.

Note: this plugin has been quickly tested, but never really used.

## Build

Launch one of the CMake project configuration tool and select your build directory. Set the CMAKE_MODULE_PATH variable to the location of the FindTULIP.cmake file (should be &lt;tulip_install_dir&gt;/share/tulip).

More informations [here](http://tulip.labri.fr/TulipDrupal/?q=node/1481).

## Use

This plugin required only one parameter, _capacity_, which must be the `IntegerProperty` containing the weights of the graph.

You will end up with one root graph with two subgraphs, being respectively the original graph and the Gomory-Hu tree.

## LICENSE

This program is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.
