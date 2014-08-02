<?php

class Site {

    public $sections = array();
    public $prev = array();
    public $next = array();
    public $last = null;
    public $toc = null;
    public $style = null;
  
    function addSection($title, $file) {
        $this->sections[$title] = $file;
        $this->prev[$title] = $this->last;
        $this->next[$title] = null;
        if ($this->last) {
             $this->next[$this->last] = $title;
        }
        $this->last   = $title; 
    }

    function top($title) {
       $s = array();
       $s[] = '<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"'.
              '"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">';
       $s[] = '<html xmlns="http://www.w3.org/1999/xhtml">';
       $s[] = '<head>';
       $s[] = '<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1" />';
       if ($this->style) {
           $s[] = '<link type="text/css" media="screen" '
                . 'rel="stylesheet" href="' . $this->style .'"/>';
       }
       $s[] = "<title>$title</title>";
       $s[] = '</head>';
       $s[] = '<body>';
       $s[] = '<h1 align="center">PSCF Manual</h2>';
       return join("\n",$s);
    }
    
    function bottom() {
       $s = array();
       $s[] = '</body>';
       $s[] = '</html>';
       return join("\n",$s);
    } 
    
    function sectionURL($title) {
        return 'sec.php?section=' . $title . '&';
    }

    function links($title, $one_file=0) {
        $s[] =  '<table width="100%">';
        $s[] =  '<tr>';
        $left = '';
        if ($this->prev[$title]) {
            $prev = $this->prev[$title];
            $ref  = $this->sections[$prev];
            if ($one_file) {
                $ref = '#' . $ref;
            } else {
                $ref = $this->sectionURL($prev);
            }
            $left = '<a href="' . $ref . '">prev</a>';
        } 
        $right = '';
        if ($this->next[$title]) {
            $next = $this->next[$title];
            $ref = $this->sections[$next];
            if ($one_file) {
                $ref = '#' . $ref;
            } else {
                $ref = $this->sectionURL($next);
            }
            $right = '<a href="' . $ref . '">next</a>';
        }
        $s[] = '<td width = "30% align="left">' . $left . '</td>';
        if ($this->toc) {
            $s[] = '<td width = "40%" align="center"> <a href="' 
                 . $this->toc . '">Table of Contents</a></td>';
        }
        $s[] = '<td width = "30%" align=right>' . $right . '</td>';
        $s[] = '</tr>';
        $s[] = '</table>';
        return join("\n",$s);
    }

    function toc() {
        $s[] = '<ul>';
        foreach ($this->sections as $title => $file) {
            $ref = $this->sectionURL($title);
            $s[] = '<li><a href="' . $ref . '">' . "$title </a></li>";
        }
        $s[] = '</ul>';
        return join("\n",$s);
    }

    function section($title) {
        $file = "{$this->sections[$title]}.html";
        print $this->top($title);
        print $this->links($title,0);
        print "<h2>$title</h2>";
        require $file;
        print $this->links($title,0);
        print $this->bottom();
    }

}

$site = new Site(); 
$site->toc   = 'toc.php';
$site->style = 'style.css';
$site->addSection('PSCF - Introduction', 'intro');
#$site->addSection('Algorithms', 'algorithms');
#$site->addSection('Limitations', 'limitations');
$site->addSection('Installation', 'install');
$site->addSection('Usage and Files', 'usage');
$site->addSection('Input Script', 'script');
$site->addSection('Field Files', 'field');
$site->addSection('Appendix: Crystal Systems', 'crystal_systems');
$site->addSection('Appendix: Space Groups', 'space_groups');
$site->addSection('Appendix: Memory Usage', 'memory');

?>
