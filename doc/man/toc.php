<?php
require 'Site.php';
print $site->top('Table of Contents');
print '<br/><br/>';
print '<h2>Table of Contents</h2>';
print $site->toc();
?>
<br/><br/> <br/>
<a href="http://www.cems.umn.edu/research/morse/code/pscf/home.php">PSCF Home Page</a>
<?php
print $site->bottom();
?>
