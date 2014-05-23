




<!DOCTYPE html>
<html class="   ">
  <head prefix="og: http://ogp.me/ns# fb: http://ogp.me/ns/fb# object: http://ogp.me/ns/object# article: http://ogp.me/ns/article# profile: http://ogp.me/ns/profile#">
    <meta charset='utf-8'>
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
    
    
    <title>COOP/typedef/constants.f90 at master · zqhuang/COOP</title>
    <link rel="search" type="application/opensearchdescription+xml" href="/opensearch.xml" title="GitHub" />
    <link rel="fluid-icon" href="https://github.com/fluidicon.png" title="GitHub" />
    <link rel="apple-touch-icon" sizes="57x57" href="/apple-touch-icon-114.png" />
    <link rel="apple-touch-icon" sizes="114x114" href="/apple-touch-icon-114.png" />
    <link rel="apple-touch-icon" sizes="72x72" href="/apple-touch-icon-144.png" />
    <link rel="apple-touch-icon" sizes="144x144" href="/apple-touch-icon-144.png" />
    <meta property="fb:app_id" content="1401488693436528"/>

      <meta content="@github" name="twitter:site" /><meta content="summary" name="twitter:card" /><meta content="zqhuang/COOP" name="twitter:title" /><meta content="COOP - Cosmology Object Oriented Package" name="twitter:description" /><meta content="https://avatars2.githubusercontent.com/u/7293247?s=400" name="twitter:image:src" />
<meta content="GitHub" property="og:site_name" /><meta content="object" property="og:type" /><meta content="https://avatars2.githubusercontent.com/u/7293247?s=400" property="og:image" /><meta content="zqhuang/COOP" property="og:title" /><meta content="https://github.com/zqhuang/COOP" property="og:url" /><meta content="COOP - Cosmology Object Oriented Package" property="og:description" />

    <link rel="assets" href="https://assets-cdn.github.com/">
    <link rel="conduit-xhr" href="https://ghconduit.com:25035/">
    <link rel="xhr-socket" href="/_sockets" />

    <meta name="msapplication-TileImage" content="/windows-tile.png" />
    <meta name="msapplication-TileColor" content="#ffffff" />
    <meta name="selected-link" value="repo_source" data-pjax-transient />
      <meta name="google-analytics" content="UA-3769691-2">

    <meta content="collector.githubapp.com" name="octolytics-host" /><meta content="collector-cdn.github.com" name="octolytics-script-host" /><meta content="github" name="octolytics-app-id" /><meta content="8E01D656:6EC0:5425114:537F9722" name="octolytics-dimension-request_id" /><meta content="7293247" name="octolytics-actor-id" /><meta content="zqhuang" name="octolytics-actor-login" /><meta content="4d1bb100bd79655f7c7cf9d5203e72b7523200d90434a1653e145199fee1ccb6" name="octolytics-actor-hash" />
    

    
    
    <link rel="icon" type="image/x-icon" href="https://assets-cdn.github.com/favicon.ico" />

    <meta content="authenticity_token" name="csrf-param" />
<meta content="sf+mipeizF5wGZvc/+euAcs7ZurWvTUAgQCNAR4IQytdX3zcbxSa9C4tjOFb3DGYn1GJeoLA2PXH1tXcavbOkA==" name="csrf-token" />

    <link href="https://assets-cdn.github.com/assets/github-a5b1ed2a59619f83449f49f1af15d646a0a9d0fd.css" media="all" rel="stylesheet" type="text/css" />
    <link href="https://assets-cdn.github.com/assets/github2-30ee11d3701e822d678d80383171fa63a23f9577.css" media="all" rel="stylesheet" type="text/css" />
    


    <meta http-equiv="x-pjax-version" content="296d4a7783891e3920ce292d67085538">

      
  <meta name="description" content="COOP - Cosmology Object Oriented Package" />

  <meta content="7293247" name="octolytics-dimension-user_id" /><meta content="zqhuang" name="octolytics-dimension-user_login" /><meta content="20108690" name="octolytics-dimension-repository_id" /><meta content="zqhuang/COOP" name="octolytics-dimension-repository_nwo" /><meta content="true" name="octolytics-dimension-repository_public" /><meta content="false" name="octolytics-dimension-repository_is_fork" /><meta content="20108690" name="octolytics-dimension-repository_network_root_id" /><meta content="zqhuang/COOP" name="octolytics-dimension-repository_network_root_nwo" />
  <link href="https://github.com/zqhuang/COOP/commits/master.atom" rel="alternate" title="Recent Commits to COOP:master" type="application/atom+xml" />

  </head>


  <body class="logged_in  env-production macintosh vis-public page-blob">
    <a href="#start-of-content" tabindex="1" class="accessibility-aid js-skip-to-content">Skip to content</a>
    <div class="wrapper">
      
      
      
      


      <div class="header header-logged-in true">
  <div class="container clearfix">

    <a class="header-logo-invertocat" href="https://github.com/">
  <span class="mega-octicon octicon-mark-github"></span>
</a>

    
    <a href="/notifications" aria-label="You have no unread notifications" class="notification-indicator tooltipped tooltipped-s" data-hotkey="g n">
        <span class="mail-status all-read"></span>
</a>

      <div class="command-bar js-command-bar  in-repository">
          <form accept-charset="UTF-8" action="/search" class="command-bar-form" id="top_search_form" method="get">

<div class="commandbar">
  <span class="message"></span>
  <input type="text" data-hotkey="s, /" name="q" id="js-command-bar-field" placeholder="Search or type a command" tabindex="1" autocapitalize="off"
    
    data-username="zqhuang"
      data-repo="zqhuang/COOP"
      data-branch="master"
      data-sha="b3e4311788eda2f86d73af0842c011e12505db3e"
  >
  <div class="display hidden"></div>
</div>

    <input type="hidden" name="nwo" value="zqhuang/COOP" />

    <div class="select-menu js-menu-container js-select-menu search-context-select-menu">
      <span class="minibutton select-menu-button js-menu-target" role="button" aria-haspopup="true">
        <span class="js-select-button">This repository</span>
      </span>

      <div class="select-menu-modal-holder js-menu-content js-navigation-container" aria-hidden="true">
        <div class="select-menu-modal">

          <div class="select-menu-item js-navigation-item js-this-repository-navigation-item selected">
            <span class="select-menu-item-icon octicon octicon-check"></span>
            <input type="radio" class="js-search-this-repository" name="search_target" value="repository" checked="checked" />
            <div class="select-menu-item-text js-select-button-text">This repository</div>
          </div> <!-- /.select-menu-item -->

          <div class="select-menu-item js-navigation-item js-all-repositories-navigation-item">
            <span class="select-menu-item-icon octicon octicon-check"></span>
            <input type="radio" name="search_target" value="global" />
            <div class="select-menu-item-text js-select-button-text">All repositories</div>
          </div> <!-- /.select-menu-item -->

        </div>
      </div>
    </div>

  <span class="help tooltipped tooltipped-s" aria-label="Show command bar help">
    <span class="octicon octicon-question"></span>
  </span>


  <input type="hidden" name="ref" value="cmdform">

</form>
        <ul class="top-nav">
          <li class="explore"><a href="/explore">Explore</a></li>
            <li><a href="https://gist.github.com">Gist</a></li>
            <li><a href="/blog">Blog</a></li>
          <li><a href="https://help.github.com">Help</a></li>
        </ul>
      </div>

    


  <ul id="user-links">
    <li>
      <a href="/zqhuang" class="name">
        <img alt="zqhuang" class=" js-avatar" data-user="7293247" height="20" src="https://avatars1.githubusercontent.com/u/7293247?s=140" width="20" /> zqhuang
      </a>
    </li>

    <li class="new-menu dropdown-toggle js-menu-container">
      <a href="#" class="js-menu-target tooltipped tooltipped-s" aria-label="Create new...">
        <span class="octicon octicon-plus"></span>
        <span class="dropdown-arrow"></span>
      </a>

      <div class="new-menu-content js-menu-content">
      </div>
    </li>

    <li>
      <a href="/settings/profile" id="account_settings"
        class="tooltipped tooltipped-s"
        aria-label="Account settings ">
        <span class="octicon octicon-tools"></span>
      </a>
    </li>
    <li>
      <form class="logout-form" action="/logout" method="post">
        <button class="sign-out-button tooltipped tooltipped-s" aria-label="Sign out">
          <span class="octicon octicon-sign-out"></span>
        </button>
      </form>
    </li>

  </ul>

<div class="js-new-dropdown-contents hidden">
  

<ul class="dropdown-menu">
  <li>
    <a href="/new"><span class="octicon octicon-repo"></span> New repository</a>
  </li>
  <li>
    <a href="/organizations/new"><span class="octicon octicon-organization"></span> New organization</a>
  </li>


    <li class="section-title">
      <span title="zqhuang/COOP">This repository</span>
    </li>
      <li>
        <a href="/zqhuang/COOP/issues/new"><span class="octicon octicon-issue-opened"></span> New issue</a>
      </li>
      <li>
        <a href="/zqhuang/COOP/settings/collaboration"><span class="octicon octicon-person"></span> New collaborator</a>
      </li>
</ul>

</div>


    
  </div>
</div>

      

        



      <div id="start-of-content" class="accessibility-aid"></div>
          <div class="site" itemscope itemtype="http://schema.org/WebPage">
    <div id="js-flash-container">
      
    </div>
    <div class="pagehead repohead instapaper_ignore readability-menu">
      <div class="container">
        

<ul class="pagehead-actions">

    <li class="subscription">
      <form accept-charset="UTF-8" action="/notifications/subscribe" class="js-social-container" data-autosubmit="true" data-remote="true" method="post"><div style="margin:0;padding:0;display:inline"><input name="authenticity_token" type="hidden" value="q1/h05PEzDvLbV3XCcIFQaG1+qgcPKvtKaT4HAjRTPkkeTnLKqgrMgPCbny0W7TWFRJV4/4Ap7WaHKtHZVvcsw==" /></div>  <input id="repository_id" name="repository_id" type="hidden" value="20108690" />

    <div class="select-menu js-menu-container js-select-menu">
      <a class="social-count js-social-count" href="/zqhuang/COOP/watchers">
        1
      </a>
      <span class="minibutton select-menu-button with-count js-menu-target" role="button" tabindex="0" aria-haspopup="true">
        <span class="js-select-button">
          <span class="octicon octicon-eye"></span>
          Unwatch
        </span>
      </span>

      <div class="select-menu-modal-holder">
        <div class="select-menu-modal subscription-menu-modal js-menu-content" aria-hidden="true">
          <div class="select-menu-header">
            <span class="select-menu-title">Notification status</span>
            <span class="octicon octicon-x js-menu-close"></span>
          </div> <!-- /.select-menu-header -->

          <div class="select-menu-list js-navigation-container" role="menu">

            <div class="select-menu-item js-navigation-item " role="menuitem" tabindex="0">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <div class="select-menu-item-text">
                <input id="do_included" name="do" type="radio" value="included" />
                <h4>Not watching</h4>
                <span class="description">You only receive notifications for conversations in which you participate or are @mentioned.</span>
                <span class="js-select-button-text hidden-select-button-text">
                  <span class="octicon octicon-eye"></span>
                  Watch
                </span>
              </div>
            </div> <!-- /.select-menu-item -->

            <div class="select-menu-item js-navigation-item selected" role="menuitem" tabindex="0">
              <span class="select-menu-item-icon octicon octicon octicon-check"></span>
              <div class="select-menu-item-text">
                <input checked="checked" id="do_subscribed" name="do" type="radio" value="subscribed" />
                <h4>Watching</h4>
                <span class="description">You receive notifications for all conversations in this repository.</span>
                <span class="js-select-button-text hidden-select-button-text">
                  <span class="octicon octicon-eye"></span>
                  Unwatch
                </span>
              </div>
            </div> <!-- /.select-menu-item -->

            <div class="select-menu-item js-navigation-item " role="menuitem" tabindex="0">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <div class="select-menu-item-text">
                <input id="do_ignore" name="do" type="radio" value="ignore" />
                <h4>Ignoring</h4>
                <span class="description">You do not receive any notifications for conversations in this repository.</span>
                <span class="js-select-button-text hidden-select-button-text">
                  <span class="octicon octicon-mute"></span>
                  Stop ignoring
                </span>
              </div>
            </div> <!-- /.select-menu-item -->

          </div> <!-- /.select-menu-list -->

        </div> <!-- /.select-menu-modal -->
      </div> <!-- /.select-menu-modal-holder -->
    </div> <!-- /.select-menu -->

</form>
    </li>

  <li>
  

  <div class="js-toggler-container js-social-container starring-container ">

    <form accept-charset="UTF-8" action="/zqhuang/COOP/unstar" class="js-toggler-form starred" data-remote="true" method="post"><div style="margin:0;padding:0;display:inline"><input name="authenticity_token" type="hidden" value="LJsvze/VN5R5KMV62nCfhBvi9IgrnCa2ICmLNpbI61C1ICPNTmDsceOGKqQVt3geVluJmtY8m7nuWccyaU7wYQ==" /></div>
      <button
        class="minibutton with-count js-toggler-target star-button"
        aria-label="Unstar this repository" title="Unstar zqhuang/COOP">
        <span class="octicon octicon-star"></span><span class="text">Unstar</span>
      </button>
        <a class="social-count js-social-count" href="/zqhuang/COOP/stargazers">
          0
        </a>
</form>
    <form accept-charset="UTF-8" action="/zqhuang/COOP/star" class="js-toggler-form unstarred" data-remote="true" method="post"><div style="margin:0;padding:0;display:inline"><input name="authenticity_token" type="hidden" value="Dfgie2ULSKo+y6TUgEHNb+pxkwVF6BN/ADEF0r4SlFkfKQcPxMV2RqZuxvQyC7lteHbZzob+t5+Gzd3iA8pCFw==" /></div>
      <button
        class="minibutton with-count js-toggler-target star-button"
        aria-label="Star this repository" title="Star zqhuang/COOP">
        <span class="octicon octicon-star"></span><span class="text">Star</span>
      </button>
        <a class="social-count js-social-count" href="/zqhuang/COOP/stargazers">
          0
        </a>
</form>  </div>

  </li>


        <li>
          <a href="/zqhuang/COOP/fork" class="minibutton with-count js-toggler-target fork-button lighter tooltipped-n" title="Fork your own copy of zqhuang/COOP to your account" aria-label="Fork your own copy of zqhuang/COOP to your account" rel="facebox nofollow">
            <span class="octicon octicon-repo-forked"></span><span class="text">Fork</span>
          </a>
          <a href="/zqhuang/COOP/network" class="social-count">0</a>
        </li>


</ul>

        <h1 itemscope itemtype="http://data-vocabulary.org/Breadcrumb" class="entry-title public">
          <span class="repo-label"><span>public</span></span>
          <span class="mega-octicon octicon-repo"></span>
          <span class="author"><a href="/zqhuang" class="url fn" itemprop="url" rel="author"><span itemprop="title">zqhuang</span></a></span><!--
       --><span class="path-divider">/</span><!--
       --><strong><a href="/zqhuang/COOP" class="js-current-repository js-repo-home-link">COOP</a></strong>

          <span class="page-context-loader">
            <img alt="" height="16" src="https://assets-cdn.github.com/images/spinners/octocat-spinner-32.gif" width="16" />
          </span>

        </h1>
      </div><!-- /.container -->
    </div><!-- /.repohead -->

    <div class="container">
      <div class="repository-with-sidebar repo-container new-discussion-timeline js-new-discussion-timeline  ">
        <div class="repository-sidebar clearfix">
            

<div class="sunken-menu vertical-right repo-nav js-repo-nav js-repository-container-pjax js-octicon-loaders">
  <div class="sunken-menu-contents">
    <ul class="sunken-menu-group">
      <li class="tooltipped tooltipped-w" aria-label="Code">
        <a href="/zqhuang/COOP" aria-label="Code" class="selected js-selected-navigation-item sunken-menu-item" data-hotkey="g c" data-pjax="true" data-selected-links="repo_source repo_downloads repo_commits repo_releases repo_tags repo_branches /zqhuang/COOP">
          <span class="octicon octicon-code"></span> <span class="full-word">Code</span>
          <img alt="" class="mini-loader" height="16" src="https://assets-cdn.github.com/images/spinners/octocat-spinner-32.gif" width="16" />
</a>      </li>

        <li class="tooltipped tooltipped-w" aria-label="Issues">
          <a href="/zqhuang/COOP/issues" aria-label="Issues" class="js-selected-navigation-item sunken-menu-item js-disable-pjax" data-hotkey="g i" data-selected-links="repo_issues /zqhuang/COOP/issues">
            <span class="octicon octicon-issue-opened"></span> <span class="full-word">Issues</span>
            <span class='counter'>0</span>
            <img alt="" class="mini-loader" height="16" src="https://assets-cdn.github.com/images/spinners/octocat-spinner-32.gif" width="16" />
</a>        </li>

      <li class="tooltipped tooltipped-w" aria-label="Pull Requests">
        <a href="/zqhuang/COOP/pulls" aria-label="Pull Requests" class="js-selected-navigation-item sunken-menu-item js-disable-pjax" data-hotkey="g p" data-selected-links="repo_pulls /zqhuang/COOP/pulls">
            <span class="octicon octicon-git-pull-request"></span> <span class="full-word">Pull Requests</span>
            <span class='counter'>0</span>
            <img alt="" class="mini-loader" height="16" src="https://assets-cdn.github.com/images/spinners/octocat-spinner-32.gif" width="16" />
</a>      </li>


        <li class="tooltipped tooltipped-w" aria-label="Wiki">
          <a href="/zqhuang/COOP/wiki" aria-label="Wiki" class="js-selected-navigation-item sunken-menu-item js-disable-pjax" data-hotkey="g w" data-selected-links="repo_wiki /zqhuang/COOP/wiki">
            <span class="octicon octicon-book"></span> <span class="full-word">Wiki</span>
            <img alt="" class="mini-loader" height="16" src="https://assets-cdn.github.com/images/spinners/octocat-spinner-32.gif" width="16" />
</a>        </li>
    </ul>
    <div class="sunken-menu-separator"></div>
    <ul class="sunken-menu-group">

      <li class="tooltipped tooltipped-w" aria-label="Pulse">
        <a href="/zqhuang/COOP/pulse" aria-label="Pulse" class="js-selected-navigation-item sunken-menu-item" data-pjax="true" data-selected-links="pulse /zqhuang/COOP/pulse">
          <span class="octicon octicon-pulse"></span> <span class="full-word">Pulse</span>
          <img alt="" class="mini-loader" height="16" src="https://assets-cdn.github.com/images/spinners/octocat-spinner-32.gif" width="16" />
</a>      </li>

      <li class="tooltipped tooltipped-w" aria-label="Graphs">
        <a href="/zqhuang/COOP/graphs" aria-label="Graphs" class="js-selected-navigation-item sunken-menu-item" data-pjax="true" data-selected-links="repo_graphs repo_contributors /zqhuang/COOP/graphs">
          <span class="octicon octicon-graph"></span> <span class="full-word">Graphs</span>
          <img alt="" class="mini-loader" height="16" src="https://assets-cdn.github.com/images/spinners/octocat-spinner-32.gif" width="16" />
</a>      </li>

      <li class="tooltipped tooltipped-w" aria-label="Network">
        <a href="/zqhuang/COOP/network" aria-label="Network" class="js-selected-navigation-item sunken-menu-item js-disable-pjax" data-selected-links="repo_network /zqhuang/COOP/network">
          <span class="octicon octicon-repo-forked"></span> <span class="full-word">Network</span>
          <img alt="" class="mini-loader" height="16" src="https://assets-cdn.github.com/images/spinners/octocat-spinner-32.gif" width="16" />
</a>      </li>
    </ul>


      <div class="sunken-menu-separator"></div>
      <ul class="sunken-menu-group">
        <li class="tooltipped tooltipped-w" aria-label="Settings">
          <a href="/zqhuang/COOP/settings" aria-label="Settings" class="js-selected-navigation-item sunken-menu-item" data-pjax="true" data-selected-links="repo_settings /zqhuang/COOP/settings">
            <span class="octicon octicon-tools"></span> <span class="full-word">Settings</span>
            <img alt="" class="mini-loader" height="16" src="https://assets-cdn.github.com/images/spinners/octocat-spinner-32.gif" width="16" />
</a>        </li>
      </ul>
  </div>
</div>

              <div class="only-with-full-nav">
                

  

<div class="clone-url open"
  data-protocol-type="http"
  data-url="/users/set_protocol?protocol_selector=http&amp;protocol_type=push">
  <h3><strong>HTTPS</strong> clone URL</h3>
  <div class="clone-url-box">
    <input type="text" class="clone js-url-field"
           value="https://github.com/zqhuang/COOP.git" readonly="readonly">
    <span class="url-box-clippy">
    <button aria-label="copy to clipboard" class="js-zeroclipboard minibutton zeroclipboard-button" data-clipboard-text="https://github.com/zqhuang/COOP.git" data-copied-hint="copied!" type="button"><span class="octicon octicon-clippy"></span></button>
    </span>
  </div>
</div>

  

<div class="clone-url "
  data-protocol-type="ssh"
  data-url="/users/set_protocol?protocol_selector=ssh&amp;protocol_type=push">
  <h3><strong>SSH</strong> clone URL</h3>
  <div class="clone-url-box">
    <input type="text" class="clone js-url-field"
           value="git@github.com:zqhuang/COOP.git" readonly="readonly">
    <span class="url-box-clippy">
    <button aria-label="copy to clipboard" class="js-zeroclipboard minibutton zeroclipboard-button" data-clipboard-text="git@github.com:zqhuang/COOP.git" data-copied-hint="copied!" type="button"><span class="octicon octicon-clippy"></span></button>
    </span>
  </div>
</div>

  

<div class="clone-url "
  data-protocol-type="subversion"
  data-url="/users/set_protocol?protocol_selector=subversion&amp;protocol_type=push">
  <h3><strong>Subversion</strong> checkout URL</h3>
  <div class="clone-url-box">
    <input type="text" class="clone js-url-field"
           value="https://github.com/zqhuang/COOP" readonly="readonly">
    <span class="url-box-clippy">
    <button aria-label="copy to clipboard" class="js-zeroclipboard minibutton zeroclipboard-button" data-clipboard-text="https://github.com/zqhuang/COOP" data-copied-hint="copied!" type="button"><span class="octicon octicon-clippy"></span></button>
    </span>
  </div>
</div>


<p class="clone-options">You can clone with
      <a href="#" class="js-clone-selector" data-protocol="http">HTTPS</a>,
      <a href="#" class="js-clone-selector" data-protocol="ssh">SSH</a>,
      or <a href="#" class="js-clone-selector" data-protocol="subversion">Subversion</a>.
  <span class="help tooltipped tooltipped-n" aria-label="Get help on which URL is right for you.">
    <a href="https://help.github.com/articles/which-remote-url-should-i-use">
    <span class="octicon octicon-question"></span>
    </a>
  </span>
</p>

  <a href="http://mac.github.com" data-url="github-mac://openRepo/https://github.com/zqhuang/COOP" class="minibutton sidebar-button js-conduit-rewrite-url" title="Save zqhuang/COOP to your computer and use it in GitHub Desktop." aria-label="Save zqhuang/COOP to your computer and use it in GitHub Desktop.">
    <span class="octicon octicon-device-desktop"></span>
    Clone in Desktop
  </a>


                <a href="/zqhuang/COOP/archive/master.zip"
                   class="minibutton sidebar-button"
                   aria-label="Download zqhuang/COOP as a zip file"
                   title="Download zqhuang/COOP as a zip file"
                   rel="nofollow">
                  <span class="octicon octicon-cloud-download"></span>
                  Download ZIP
                </a>
              </div>
        </div><!-- /.repository-sidebar -->

        <div id="js-repo-pjax-container" class="repository-content context-loader-container" data-pjax-container>
          


<a href="/zqhuang/COOP/blob/bac85cfd65bf47d791e9b6c068b13314c55af22e/typedef/constants.f90" class="hidden js-permalink-shortcut" data-hotkey="y">Permalink</a>

<!-- blob contrib key: blob_contributors:v21:b3f6ccd97581938399806a1af308a29e -->

<p title="This is a placeholder element" class="js-history-link-replace hidden"></p>

<a href="/zqhuang/COOP/find/master" data-pjax data-hotkey="t" class="js-show-file-finder" style="display:none">Show File Finder</a>

<div class="file-navigation">
  

<div class="select-menu js-menu-container js-select-menu" >
  <span class="minibutton select-menu-button js-menu-target" data-hotkey="w"
    data-master-branch="master"
    data-ref="master"
    role="button" aria-label="Switch branches or tags" tabindex="0" aria-haspopup="true">
    <span class="octicon octicon-git-branch"></span>
    <i>branch:</i>
    <span class="js-select-button">master</span>
  </span>

  <div class="select-menu-modal-holder js-menu-content js-navigation-container" data-pjax aria-hidden="true">

    <div class="select-menu-modal">
      <div class="select-menu-header">
        <span class="select-menu-title">Switch branches/tags</span>
        <span class="octicon octicon-x js-menu-close"></span>
      </div> <!-- /.select-menu-header -->

      <div class="select-menu-filters">
        <div class="select-menu-text-filter">
          <input type="text" aria-label="Find or create a branch…" id="context-commitish-filter-field" class="js-filterable-field js-navigation-enable" placeholder="Find or create a branch…">
        </div>
        <div class="select-menu-tabs">
          <ul>
            <li class="select-menu-tab">
              <a href="#" data-tab-filter="branches" class="js-select-menu-tab">Branches</a>
            </li>
            <li class="select-menu-tab">
              <a href="#" data-tab-filter="tags" class="js-select-menu-tab">Tags</a>
            </li>
          </ul>
        </div><!-- /.select-menu-tabs -->
      </div><!-- /.select-menu-filters -->

      <div class="select-menu-list select-menu-tab-bucket js-select-menu-tab-bucket" data-tab-filter="branches">

        <div data-filterable-for="context-commitish-filter-field" data-filterable-type="substring">


            <div class="select-menu-item js-navigation-item selected">
              <span class="select-menu-item-icon octicon octicon-check"></span>
              <a href="/zqhuang/COOP/blob/master/typedef/constants.f90"
                 data-name="master"
                 data-skip-pjax="true"
                 rel="nofollow"
                 class="js-navigation-open select-menu-item-text js-select-button-text css-truncate-target"
                 title="master">master</a>
            </div> <!-- /.select-menu-item -->
        </div>

          <form accept-charset="UTF-8" action="/zqhuang/COOP/branches" class="js-create-branch select-menu-item select-menu-new-item-form js-navigation-item js-new-item-form" method="post"><div style="margin:0;padding:0;display:inline"><input name="authenticity_token" type="hidden" value="GLDQdwtfSpLBpsvN5eJxsWJvKXg3KqSxLoN5qUFcklR+i4PMIFvKmXfan7jXiSed6lIKFlyQ9Dwf6piOrDutLw==" /></div>
            <span class="octicon octicon-git-branch select-menu-item-icon"></span>
            <div class="select-menu-item-text">
              <h4>Create branch: <span class="js-new-item-name"></span></h4>
              <span class="description">from ‘master’</span>
            </div>
            <input type="hidden" name="name" id="name" class="js-new-item-value">
            <input type="hidden" name="branch" id="branch" value="master" />
            <input type="hidden" name="path" id="path" value="typedef/constants.f90" />
          </form> <!-- /.select-menu-item -->

      </div> <!-- /.select-menu-list -->

      <div class="select-menu-list select-menu-tab-bucket js-select-menu-tab-bucket" data-tab-filter="tags">
        <div data-filterable-for="context-commitish-filter-field" data-filterable-type="substring">


        </div>

        <div class="select-menu-no-results">Nothing to show</div>
      </div> <!-- /.select-menu-list -->

    </div> <!-- /.select-menu-modal -->
  </div> <!-- /.select-menu-modal-holder -->
</div> <!-- /.select-menu -->

  <div class="breadcrumb">
    <span class='repo-root js-repo-root'><span itemscope="" itemtype="http://data-vocabulary.org/Breadcrumb"><a href="/zqhuang/COOP" data-branch="master" data-direction="back" data-pjax="true" itemscope="url"><span itemprop="title">COOP</span></a></span></span><span class="separator"> / </span><span itemscope="" itemtype="http://data-vocabulary.org/Breadcrumb"><a href="/zqhuang/COOP/tree/master/typedef" data-branch="master" data-direction="back" data-pjax="true" itemscope="url"><span itemprop="title">typedef</span></a></span><span class="separator"> / </span><strong class="final-path">constants.f90</strong> <button aria-label="copy to clipboard" class="js-zeroclipboard minibutton zeroclipboard-button" data-clipboard-text="typedef/constants.f90" data-copied-hint="copied!" type="button"><span class="octicon octicon-clippy"></span></button>
  </div>
</div>


  <div class="commit file-history-tease">
      <img alt="zqhuang" class="main-avatar js-avatar" data-user="7293247" height="24" src="https://avatars1.githubusercontent.com/u/7293247?s=140" width="24" />
      <span class="author"><a href="/zqhuang" rel="author">zqhuang</a></span>
      <time datetime="2014-05-23T14:28:13-04:00" is="relative-time" title-format="%Y-%m-%d %H:%M:%S %z" title="2014-05-23 14:28:13 -0400">May 23, 2014</time>
      <div class="commit-title">
          <a href="/zqhuang/COOP/commit/24488b94a32f35d9ee2dd11c7c84c85a392cc5c0" class="message" data-pjax="true" title="set up structure

basic structure">set up structure</a>
      </div>

    <div class="participation">
      <p class="quickstat"><a href="#blob_contributors_box" rel="facebox"><strong>1</strong>  contributor</a></p>
      
    </div>
    <div id="blob_contributors_box" style="display:none">
      <h2 class="facebox-header">Users who have contributed to this file</h2>
      <ul class="facebox-user-list">
          <li class="facebox-user-list-item">
            <img alt="zqhuang" class=" js-avatar" data-user="7293247" height="24" src="https://avatars1.githubusercontent.com/u/7293247?s=140" width="24" />
            <a href="/zqhuang">zqhuang</a>
          </li>
      </ul>
    </div>
  </div>

<div class="file-box">
  <div class="file">
    <div class="meta clearfix">
      <div class="info file-name">
        <span class="icon"><b class="octicon octicon-file-text"></b></span>
        <span class="mode" title="File Mode">file</span>
        <span class="meta-divider"></span>
          <span>128 lines (108 sloc)</span>
          <span class="meta-divider"></span>
        <span>7.534 kb</span>
      </div>
      <div class="actions">
        <div class="button-group">
            <a class="minibutton tooltipped tooltipped-w js-conduit-openfile-check"
               href="http://mac.github.com"
               data-url="github-mac://openRepo/https://github.com/zqhuang/COOP?branch=master&amp;filepath=typedef%2Fconstants.f90"
               aria-label="Open this file in GitHub for Mac"
               data-failed-title="Your version of GitHub for Mac is too old to open this file. Try checking for updates.">
                <span class="octicon octicon-device-desktop"></span> Open
            </a>
                <a class="minibutton js-update-url-with-hash"
                   href="/zqhuang/COOP/edit/master/typedef/constants.f90"
                   data-method="post" rel="nofollow" data-hotkey="e">Edit</a>
          <a href="/zqhuang/COOP/raw/master/typedef/constants.f90" class="button minibutton " id="raw-url">Raw</a>
            <a href="/zqhuang/COOP/blame/master/typedef/constants.f90" class="button minibutton js-update-url-with-hash">Blame</a>
          <a href="/zqhuang/COOP/commits/master/typedef/constants.f90" class="button minibutton " rel="nofollow">History</a>
        </div><!-- /.button-group -->

            <a class="minibutton danger empty-icon"
               href="/zqhuang/COOP/delete/master/typedef/constants.f90"
               data-method="post" data-test-id="delete-blob-file" rel="nofollow">

          Delete
        </a>
      </div><!-- /.actions -->
    </div>
        <div class="blob-wrapper data type-fortran js-blob-data">
        <table class="file-code file-diff tab-size-8">
          <tr class="file-code-line">
            <td class="blob-line-nums">
              <span id="L1" rel="#L1">1</span>
<span id="L2" rel="#L2">2</span>
<span id="L3" rel="#L3">3</span>
<span id="L4" rel="#L4">4</span>
<span id="L5" rel="#L5">5</span>
<span id="L6" rel="#L6">6</span>
<span id="L7" rel="#L7">7</span>
<span id="L8" rel="#L8">8</span>
<span id="L9" rel="#L9">9</span>
<span id="L10" rel="#L10">10</span>
<span id="L11" rel="#L11">11</span>
<span id="L12" rel="#L12">12</span>
<span id="L13" rel="#L13">13</span>
<span id="L14" rel="#L14">14</span>
<span id="L15" rel="#L15">15</span>
<span id="L16" rel="#L16">16</span>
<span id="L17" rel="#L17">17</span>
<span id="L18" rel="#L18">18</span>
<span id="L19" rel="#L19">19</span>
<span id="L20" rel="#L20">20</span>
<span id="L21" rel="#L21">21</span>
<span id="L22" rel="#L22">22</span>
<span id="L23" rel="#L23">23</span>
<span id="L24" rel="#L24">24</span>
<span id="L25" rel="#L25">25</span>
<span id="L26" rel="#L26">26</span>
<span id="L27" rel="#L27">27</span>
<span id="L28" rel="#L28">28</span>
<span id="L29" rel="#L29">29</span>
<span id="L30" rel="#L30">30</span>
<span id="L31" rel="#L31">31</span>
<span id="L32" rel="#L32">32</span>
<span id="L33" rel="#L33">33</span>
<span id="L34" rel="#L34">34</span>
<span id="L35" rel="#L35">35</span>
<span id="L36" rel="#L36">36</span>
<span id="L37" rel="#L37">37</span>
<span id="L38" rel="#L38">38</span>
<span id="L39" rel="#L39">39</span>
<span id="L40" rel="#L40">40</span>
<span id="L41" rel="#L41">41</span>
<span id="L42" rel="#L42">42</span>
<span id="L43" rel="#L43">43</span>
<span id="L44" rel="#L44">44</span>
<span id="L45" rel="#L45">45</span>
<span id="L46" rel="#L46">46</span>
<span id="L47" rel="#L47">47</span>
<span id="L48" rel="#L48">48</span>
<span id="L49" rel="#L49">49</span>
<span id="L50" rel="#L50">50</span>
<span id="L51" rel="#L51">51</span>
<span id="L52" rel="#L52">52</span>
<span id="L53" rel="#L53">53</span>
<span id="L54" rel="#L54">54</span>
<span id="L55" rel="#L55">55</span>
<span id="L56" rel="#L56">56</span>
<span id="L57" rel="#L57">57</span>
<span id="L58" rel="#L58">58</span>
<span id="L59" rel="#L59">59</span>
<span id="L60" rel="#L60">60</span>
<span id="L61" rel="#L61">61</span>
<span id="L62" rel="#L62">62</span>
<span id="L63" rel="#L63">63</span>
<span id="L64" rel="#L64">64</span>
<span id="L65" rel="#L65">65</span>
<span id="L66" rel="#L66">66</span>
<span id="L67" rel="#L67">67</span>
<span id="L68" rel="#L68">68</span>
<span id="L69" rel="#L69">69</span>
<span id="L70" rel="#L70">70</span>
<span id="L71" rel="#L71">71</span>
<span id="L72" rel="#L72">72</span>
<span id="L73" rel="#L73">73</span>
<span id="L74" rel="#L74">74</span>
<span id="L75" rel="#L75">75</span>
<span id="L76" rel="#L76">76</span>
<span id="L77" rel="#L77">77</span>
<span id="L78" rel="#L78">78</span>
<span id="L79" rel="#L79">79</span>
<span id="L80" rel="#L80">80</span>
<span id="L81" rel="#L81">81</span>
<span id="L82" rel="#L82">82</span>
<span id="L83" rel="#L83">83</span>
<span id="L84" rel="#L84">84</span>
<span id="L85" rel="#L85">85</span>
<span id="L86" rel="#L86">86</span>
<span id="L87" rel="#L87">87</span>
<span id="L88" rel="#L88">88</span>
<span id="L89" rel="#L89">89</span>
<span id="L90" rel="#L90">90</span>
<span id="L91" rel="#L91">91</span>
<span id="L92" rel="#L92">92</span>
<span id="L93" rel="#L93">93</span>
<span id="L94" rel="#L94">94</span>
<span id="L95" rel="#L95">95</span>
<span id="L96" rel="#L96">96</span>
<span id="L97" rel="#L97">97</span>
<span id="L98" rel="#L98">98</span>
<span id="L99" rel="#L99">99</span>
<span id="L100" rel="#L100">100</span>
<span id="L101" rel="#L101">101</span>
<span id="L102" rel="#L102">102</span>
<span id="L103" rel="#L103">103</span>
<span id="L104" rel="#L104">104</span>
<span id="L105" rel="#L105">105</span>
<span id="L106" rel="#L106">106</span>
<span id="L107" rel="#L107">107</span>
<span id="L108" rel="#L108">108</span>
<span id="L109" rel="#L109">109</span>
<span id="L110" rel="#L110">110</span>
<span id="L111" rel="#L111">111</span>
<span id="L112" rel="#L112">112</span>
<span id="L113" rel="#L113">113</span>
<span id="L114" rel="#L114">114</span>
<span id="L115" rel="#L115">115</span>
<span id="L116" rel="#L116">116</span>
<span id="L117" rel="#L117">117</span>
<span id="L118" rel="#L118">118</span>
<span id="L119" rel="#L119">119</span>
<span id="L120" rel="#L120">120</span>
<span id="L121" rel="#L121">121</span>
<span id="L122" rel="#L122">122</span>
<span id="L123" rel="#L123">123</span>
<span id="L124" rel="#L124">124</span>
<span id="L125" rel="#L125">125</span>
<span id="L126" rel="#L126">126</span>
<span id="L127" rel="#L127">127</span>

            </td>
            <td class="blob-line-code"><div class="code-body highlight"><pre><div class='line' id='LC1'><span class="k">module </span><span class="nv">coop_constants</span></div><div class='line' id='LC2'>&nbsp;&nbsp;<span class="k">implicit none</span></div><div class='line' id='LC3'><br/></div><div class='line' id='LC4'><span class="err">#</span><span class="k">include</span> <span class="s2">&quot;constants.h&quot;</span></div><div class='line' id='LC5'>&nbsp;&nbsp;<span class="kt">integer</span><span class="p">,</span><span class="k">parameter</span><span class="kd">::</span><span class="nv">coop_real_length</span> <span class="o">=</span> <span class="mi">8</span></div><div class='line' id='LC6'>&nbsp;&nbsp;<span class="kt">integer</span><span class="p">,</span><span class="k">parameter</span><span class="kd">::</span><span class="nv">coop_short_string_length</span> <span class="o">=</span> <span class="mi">32</span></div><div class='line' id='LC7'>&nbsp;&nbsp;<span class="kt">integer</span><span class="p">,</span><span class="k">parameter</span><span class="kd">::</span><span class="nv">coop_string_length</span> <span class="o">=</span> <span class="mi">256</span></div><div class='line' id='LC8'>&nbsp;&nbsp;<span class="kt">integer</span><span class="p">,</span><span class="k">parameter</span><span class="kd">::</span><span class="nv">coop_long_string_length</span> <span class="o">=</span> <span class="mi">8192</span></div><div class='line' id='LC9'><br/></div><div class='line' id='LC10'>&nbsp;&nbsp;<span class="kt">Character</span><span class="p">,</span><span class="k">Parameter</span><span class="kd">::</span><span class="nv">coop_backslash</span> <span class="o">=</span> <span class="nb">Char</span><span class="p">(</span><span class="mi">92</span><span class="p">)</span></div><div class='line' id='LC11'>&nbsp;&nbsp;<span class="kt">character</span><span class="p">,</span><span class="k">parameter</span><span class="kd">::</span><span class="nv">coop_slash</span> <span class="o">=</span> <span class="nb">Char</span><span class="p">(</span><span class="mi">47</span><span class="p">)</span></div><div class='line' id='LC12'>&nbsp;&nbsp;<span class="kt">Character</span><span class="p">,</span><span class="k">Parameter</span><span class="kd">::</span><span class="nv">coop_backspace</span> <span class="o">=</span> <span class="nb">Char</span><span class="p">(</span><span class="mi">8</span><span class="p">)</span></div><div class='line' id='LC13'>&nbsp;&nbsp;<span class="kt">Character</span><span class="p">,</span><span class="k">Parameter</span><span class="kd">::</span><span class="nv">coop_tab</span> <span class="o">=</span> <span class="nb">Char</span><span class="p">(</span><span class="mi">9</span><span class="p">)</span></div><div class='line' id='LC14'>&nbsp;&nbsp;<span class="kt">Character</span><span class="p">,</span><span class="k">Parameter</span><span class="kd">::</span><span class="nv">coop_newline</span> <span class="o">=</span> <span class="nb">Char</span><span class="p">(</span><span class="mi">10</span><span class="p">)</span></div><div class='line' id='LC15'>&nbsp;&nbsp;<span class="kt">Character</span><span class="p">,</span><span class="k">Parameter</span><span class="kd">::</span><span class="nv">coop_vertical_tab</span> <span class="o">=</span> <span class="nb">Char</span><span class="p">(</span><span class="mi">11</span><span class="p">)</span></div><div class='line' id='LC16'>&nbsp;&nbsp;<span class="kt">Character</span><span class="p">,</span><span class="k">Parameter</span><span class="kd">::</span><span class="nv">coop_newpage</span> <span class="o">=</span> <span class="nb">Char</span><span class="p">(</span><span class="mi">12</span><span class="p">)</span></div><div class='line' id='LC17'>&nbsp;&nbsp;<span class="kt">Character</span><span class="p">,</span><span class="k">Parameter</span><span class="kd">::</span><span class="nv">coop_carriage_return</span> <span class="o">=</span> <span class="nb">Char</span><span class="p">(</span><span class="mi">13</span><span class="p">)</span></div><div class='line' id='LC18'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_pi</span> <span class="o">=</span> <span class="mf">3.14159265358979323846264338327950288</span><span class="nv">d0</span></div><div class='line' id='LC19'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_pio2</span> <span class="o">=</span> <span class="nv">coop_pi</span><span class="o">/</span><span class="mf">2.</span><span class="err">_</span><span class="nv">dl</span></div><div class='line' id='LC20'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_pio4</span> <span class="o">=</span> <span class="nv">coop_pi</span><span class="o">/</span><span class="mf">4.</span><span class="err">_</span><span class="nv">dl</span></div><div class='line' id='LC21'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_ln2</span> <span class="o">=</span> <span class="mf">0.6931471805599453094172321</span><span class="nv">d0</span> </div><div class='line' id='LC22'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_ln10</span> <span class="o">=</span> <span class="mf">2.302585092994045684017991</span><span class="nv">d0</span></div><div class='line' id='LC23'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">PARAMETER</span><span class="kd">::</span> <span class="nv">coop_LnPi</span><span class="o">=</span><span class="mf">1.144729885849400174143427351</span><span class="nv">d0</span></div><div class='line' id='LC24'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">PARAMETER</span><span class="kd">::</span> <span class="nv">coop_LogPi</span><span class="o">=</span><span class="mf">0.4971498726941338543512682883</span><span class="nv">d0</span></div><div class='line' id='LC25'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_sqrt2</span> <span class="o">=</span> <span class="mf">1.4142135623730950488016887</span><span class="nv">d0</span></div><div class='line' id='LC26'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">PARAMETER</span><span class="kd">::</span> <span class="nv">coop_sqrt3</span> <span class="o">=</span> <span class="mf">1.73205080756887729352744634</span><span class="nv">d0</span></div><div class='line' id='LC27'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_sqrt5</span> <span class="o">=</span> <span class="mf">2.236067977499789696409174</span><span class="nv">d0</span></div><div class='line' id='LC28'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_sqrt6</span> <span class="o">=</span> <span class="nv">coop_sqrt2</span><span class="o">*</span><span class="nv">coop_sqrt3</span></div><div class='line' id='LC29'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_sqrt7</span> <span class="o">=</span> <span class="mf">2.645751311064590590502</span><span class="nv">d0</span></div><div class='line' id='LC30'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_sqrt8</span> <span class="o">=</span> <span class="mf">2.</span><span class="nv">d0</span><span class="o">*</span><span class="nv">coop_sqrt2</span></div><div class='line' id='LC31'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_sqrt10</span> <span class="o">=</span> <span class="nv">coop_sqrt2</span> <span class="o">*</span> <span class="nv">coop_sqrt5</span></div><div class='line' id='LC32'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_sqrt11</span> <span class="o">=</span> <span class="mf">3.316624790355399849115</span><span class="nv">d0</span></div><div class='line' id='LC33'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_sqrt12</span> <span class="o">=</span> <span class="mf">2.</span><span class="nv">d0</span><span class="o">*</span><span class="nv">coop_sqrt3</span></div><div class='line' id='LC34'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_sqrt13</span> <span class="o">=</span> <span class="mf">3.605551275463989293119</span><span class="nv">d0</span></div><div class='line' id='LC35'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_sqrt14</span> <span class="o">=</span> <span class="nv">coop_sqrt2</span><span class="o">*</span><span class="nv">coop_sqrt7</span></div><div class='line' id='LC36'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_sqrt15</span> <span class="o">=</span> <span class="nv">coop_sqrt3</span><span class="o">*</span><span class="nv">coop_sqrt5</span></div><div class='line' id='LC37'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_8pi</span> <span class="o">=</span> <span class="nv">coop_pi</span><span class="o">*</span><span class="mf">8.</span><span class="err">_</span><span class="nv">dl</span></div><div class='line' id='LC38'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_4pi</span> <span class="o">=</span> <span class="nv">coop_pi</span><span class="o">*</span><span class="mf">4.</span><span class="err">_</span><span class="nv">dl</span></div><div class='line' id='LC39'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_4piby3</span> <span class="o">=</span> <span class="nv">coop_4pi</span><span class="o">/</span><span class="mf">3.</span><span class="err">_</span><span class="nv">dl</span></div><div class='line' id='LC40'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_3by4pi</span> <span class="o">=</span> <span class="mf">3.</span><span class="err">_</span><span class="nv">dl</span><span class="o">/</span><span class="nv">coop_4pi</span></div><div class='line' id='LC41'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_2pi</span> <span class="o">=</span> <span class="nv">coop_pi</span><span class="o">*</span><span class="mf">2.</span><span class="err">_</span><span class="nv">dl</span></div><div class='line' id='LC42'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_pi2</span> <span class="o">=</span> <span class="nv">coop_pi</span> <span class="o">**</span> <span class="mi">2</span></div><div class='line' id='LC43'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_pi4</span> <span class="o">=</span> <span class="nv">coop_pi2</span> <span class="o">**</span> <span class="mi">2</span></div><div class='line' id='LC44'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_2pi2</span> <span class="o">=</span> <span class="nv">coop_pi2</span> <span class="o">*</span> <span class="mf">2.</span><span class="err">_</span><span class="nv">dl</span></div><div class='line' id='LC45'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_8pi3</span> <span class="o">=</span> <span class="nv">coop_2pi</span> <span class="o">**</span> <span class="mi">3</span></div><div class='line' id='LC46'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_7pi4by120</span> <span class="o">=</span> <span class="nv">coop_pi4</span> <span class="o">*</span> <span class="p">(</span><span class="mf">7.</span><span class="nv">d0</span><span class="o">/</span><span class="mi">12</span><span class="mf">0.</span><span class="nv">d0</span><span class="p">)</span></div><div class='line' id='LC47'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_third</span> <span class="o">=</span> <span class="mf">1.</span><span class="err">_</span><span class="nv">dl</span><span class="o">/</span><span class="mf">3.</span><span class="err">_</span><span class="nv">dl</span></div><div class='line' id='LC48'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_two_thirds</span> <span class="o">=</span> <span class="mf">2.</span><span class="err">_</span><span class="nv">dl</span><span class="o">/</span><span class="mf">3.</span><span class="err">_</span><span class="nv">dl</span></div><div class='line' id='LC49'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_four_thirds</span> <span class="o">=</span> <span class="mf">4.</span><span class="err">_</span><span class="nv">dl</span><span class="o">/</span><span class="mf">3.</span><span class="err">_</span><span class="nv">dl</span></div><div class='line' id='LC50'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_sqrtpi</span> <span class="o">=</span> <span class="mf">1.7724538509055160272981674833411</span><span class="err">_</span><span class="nv">dl</span></div><div class='line' id='LC51'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_EulerC</span><span class="o">=</span><span class="mf">0.57721566490153286060651209</span><span class="err">_</span><span class="nv">dl</span></div><div class='line' id='LC52'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_Riemannzeta3</span> <span class="o">=</span> <span class="mf">1.2020569031595942853997</span><span class="err">_</span><span class="nv">dl</span></div><div class='line' id='LC53'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_Riemannzeta5</span>  <span class="o">=</span> <span class="mf">1.0369277551433699263313</span><span class="err">_</span><span class="nv">dl</span></div><div class='line' id='LC54'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_Riemannzeta7</span>  <span class="o">=</span> <span class="mf">1.0083492773819228268397</span><span class="err">_</span><span class="nv">dl</span></div><div class='line' id='LC55'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_Riemannzeta9</span> <span class="o">=</span> <span class="mf">1.00200839282608221441785</span><span class="err">_</span><span class="nv">dl</span></div><div class='line' id='LC56'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">PARAMETER</span><span class="kd">::</span> <span class="nv">coop_fullsky_degrees</span> <span class="o">=</span> <span class="mi">4125</span><span class="mf">2.96125</span>   <span class="c">!!4pi/degree^2</span></div><div class='line' id='LC57'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_degree</span> <span class="o">=</span> <span class="nv">coop_pi</span><span class="o">/</span><span class="mi">18</span><span class="mf">0.</span><span class="nv">d0</span></div><div class='line' id='LC58'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_arcmin</span> <span class="o">=</span> <span class="nv">coop_degree</span><span class="o">/</span><span class="mi">6</span><span class="mf">0.</span><span class="nv">d0</span></div><div class='line' id='LC59'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_sigmabyfwhm</span> <span class="o">=</span> <span class="mf">1.</span><span class="nv">d0</span><span class="o">/</span><span class="nb">sqrt</span><span class="p">(</span><span class="mf">8.</span><span class="nv">d0</span><span class="o">*</span><span class="nv">coop_ln2</span><span class="p">)</span></div><div class='line' id='LC60'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_arcsec</span> <span class="o">=</span> <span class="nv">coop_arcmin</span><span class="o">/</span><span class="mi">6</span><span class="mf">0.</span><span class="nv">d0</span></div><div class='line' id='LC61'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_chbyMpcH0</span> <span class="o">=</span> <span class="mi">299</span><span class="mf">7.92458</span><span class="nv">d0</span></div><div class='line' id='LC62'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span> <span class="kd">::</span> <span class="nv">coop_NewtonG</span> <span class="o">=</span> <span class="mf">6.67428e-11</span><span class="err">_</span><span class="nv">dl</span></div><div class='line' id='LC63'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_planck</span> <span class="o">=</span> <span class="mf">6.626068</span><span class="nv">d</span><span class="o">-</span><span class="mi">34</span></div><div class='line' id='LC64'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_hbar</span> <span class="o">=</span> <span class="nv">coop_planck</span><span class="o">/</span><span class="nv">coop_2pi</span></div><div class='line' id='LC65'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_boltzmann</span> <span class="o">=</span> <span class="mf">1.3806504</span><span class="nv">d</span><span class="o">-</span><span class="mi">23</span></div><div class='line' id='LC66'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_K2GHz</span> <span class="o">=</span> <span class="nv">coop_boltzmann</span><span class="o">/</span><span class="p">(</span><span class="nv">coop_planck</span><span class="o">*</span><span class="mf">1.</span><span class="nv">d9</span><span class="p">)</span> <span class="c">!!</span></div><div class='line' id='LC67'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_Mpc_SI</span> <span class="o">=</span> <span class="mf">3.08568025e22</span><span class="err">_</span><span class="nv">dl</span></div><div class='line' id='LC68'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_Msun_SI</span> <span class="o">=</span> <span class="mf">1.98892e30</span><span class="err">_</span><span class="nv">dl</span></div><div class='line' id='LC69'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_SpeedOfLight_SI</span> <span class="o">=</span> <span class="mi">29979245</span><span class="mf">8.</span><span class="err">_</span><span class="nv">dl</span></div><div class='line' id='LC70'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_c</span> <span class="o">=</span> <span class="nv">coop_SpeedOfLight_SI</span></div><div class='line' id='LC71'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_PlanckMass_SI</span> <span class="o">=</span> <span class="mf">2.17644e-8</span><span class="err">_</span><span class="nv">dl</span></div><div class='line' id='LC72'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_PlanckLength_SI</span> <span class="o">=</span><span class="mf">1.616252e-35</span><span class="err">_</span><span class="nv">dl</span></div><div class='line' id='LC73'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_PlanckTime_SI</span> <span class="o">=</span> <span class="mf">5.39124e-44</span><span class="err">_</span><span class="nv">dl</span></div><div class='line' id='LC74'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_PlanckEnergy_SI</span> <span class="o">=</span> <span class="nv">coop_PlanckMass_SI</span> <span class="o">*</span> <span class="nv">coop_SpeedOfLight_SI</span> <span class="o">**</span> <span class="mi">2</span></div><div class='line' id='LC75'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_PlanckTemperature_SI</span> <span class="o">=</span> <span class="nv">coop_PlanckEnergy_SI</span> <span class="o">/</span> <span class="nv">coop_boltzmann</span></div><div class='line' id='LC76'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_Yr_SI</span> <span class="o">=</span> <span class="mi">360</span><span class="mf">0.</span><span class="err">_</span><span class="nv">dl</span><span class="o">*</span><span class="mi">2</span><span class="mf">4.</span><span class="err">_</span><span class="nv">dl</span><span class="o">*</span><span class="mi">36</span><span class="mf">5.2422</span><span class="err">_</span><span class="nv">dl</span> <span class="c">!!time</span></div><div class='line' id='LC77'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_Gyr_SI</span> <span class="o">=</span> <span class="nv">coop_Yr_SI</span> <span class="o">*</span> <span class="mf">1.e9</span><span class="err">_</span><span class="nv">dl</span> <span class="c">!!time</span></div><div class='line' id='LC78'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_Lyr_SI</span> <span class="o">=</span> <span class="mf">9.4605284e15</span><span class="err">_</span><span class="nv">dl</span> <span class="c">!!length</span></div><div class='line' id='LC79'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_eV_SI</span> <span class="o">=</span> <span class="mf">1.60217648740e-19</span><span class="err">_</span><span class="nv">dl</span></div><div class='line' id='LC80'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_eVmass_SI</span> <span class="o">=</span> <span class="nv">coop_eV_SI</span> <span class="o">/</span> <span class="nv">coop_SpeedOfLight_SI</span> <span class="o">**</span> <span class="mi">2</span></div><div class='line' id='LC81'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_eVlength_SI</span> <span class="o">=</span> <span class="nv">coop_SpeedOfLight_SI</span><span class="o">/</span> <span class="p">(</span><span class="nv">coop_eV_SI</span> <span class="o">/</span><span class="nv">coop_hbar</span><span class="p">)</span></div><div class='line' id='LC82'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_MeV_SI</span> <span class="o">=</span> <span class="nv">coop_eV_SI</span> <span class="o">*</span> <span class="mf">1.e6</span><span class="err">_</span><span class="nv">dl</span></div><div class='line' id='LC83'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_MeVmass_SI</span> <span class="o">=</span> <span class="nv">coop_MeV_SI</span> <span class="o">/</span> <span class="nv">coop_SpeedOfLight_SI</span> <span class="o">**</span> <span class="mi">2</span></div><div class='line' id='LC84'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_GeV_SI</span> <span class="o">=</span> <span class="nv">coop_eV_SI</span> <span class="o">*</span> <span class="mf">1.e9</span><span class="err">_</span><span class="nv">dl</span></div><div class='line' id='LC85'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_GeVmass_SI</span> <span class="o">=</span> <span class="nv">coop_GeV_SI</span> <span class="o">/</span> <span class="nv">coop_SpeedOfLight_SI</span> <span class="o">**</span> <span class="mi">2</span></div><div class='line' id='LC86'><br/></div><div class='line' id='LC87'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_fine_structure</span> <span class="o">=</span> <span class="mf">1.</span><span class="err">_</span><span class="nv">dl</span><span class="o">/</span><span class="mi">13</span><span class="mf">7.035</span><span class="err">_</span><span class="nv">dl</span></div><div class='line' id='LC88'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_atomic_mass_unit_SI</span> <span class="o">=</span> <span class="mf">1.66053878283e-27</span><span class="err">_</span><span class="nv">dl</span></div><div class='line' id='LC89'><br/></div><div class='line' id='LC90'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_electron_mass_SI</span> <span class="o">=</span> <span class="nv">coop_atomic_mass_unit_SI</span> <span class="o">/</span> <span class="mi">182</span><span class="mf">2.8884845</span><span class="err">_</span><span class="nv">dl</span> <span class="c">!!0.51099892811 * coop_MeVmass_SI</span></div><div class='line' id='LC91'><br/></div><div class='line' id='LC92'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_proton_mass_SI</span> <span class="o">=</span> <span class="mf">1.00727646681290</span> <span class="o">*</span> <span class="nv">coop_atomic_mass_unit_SI</span></div><div class='line' id='LC93'><br/></div><div class='line' id='LC94'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_E_hydrogen</span> <span class="o">=</span> <span class="mi">1</span><span class="mf">3.605698</span> <span class="o">*</span> <span class="nv">coop_eV_SI</span> </div><div class='line' id='LC95'><br/></div><div class='line' id='LC96'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_Hydrogen_mass_SI</span> <span class="o">=</span> <span class="nv">coop_electron_mass_SI</span> <span class="o">+</span> <span class="nv">coop_proton_mass_SI</span> <span class="o">-</span> <span class="nv">coop_E_Hydrogen</span><span class="o">/</span><span class="nv">coop_SpeedOfLight_SI</span><span class="o">**</span><span class="mi">2</span></div><div class='line' id='LC97'><br/></div><div class='line' id='LC98'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_massratio_He_H</span> <span class="o">=</span>  <span class="mf">3.9715</span><span class="err">_</span><span class="nv">dl</span></div><div class='line' id='LC99'><br/></div><div class='line' id='LC100'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_rhocritbyh2_SI</span> <span class="o">=</span> <span class="mf">3.</span><span class="err">_</span><span class="nv">dl</span><span class="o">*</span><span class="p">(</span><span class="mf">1.e5</span><span class="err">_</span><span class="nv">dl</span><span class="o">/</span><span class="nv">coop_SpeedOfLight_SI</span><span class="p">)</span><span class="o">**</span><span class="mi">2</span><span class="o">*</span><span class="p">(</span><span class="nv">coop_PlanckMass_SI</span><span class="o">/</span><span class="nv">coop_Mpc_SI</span><span class="o">**</span><span class="mi">2</span><span class="o">/</span><span class="nv">coop_PlanckLength_SI</span><span class="p">)</span><span class="o">/</span><span class="nv">coop_8pi</span> <span class="c">!!1.878e-26 !!mass/volume, note that this is NOT energy per volume</span></div><div class='line' id='LC101'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_hbyH0_SI</span> <span class="o">=</span> <span class="nv">coop_Mpc_SI</span> <span class="o">/</span> <span class="mf">1.e5</span><span class="err">_</span><span class="nv">dl</span> <span class="c">!!Cosmic age unit</span></div><div class='line' id='LC102'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_hbyH0_Gyr</span> <span class="o">=</span> <span class="nv">coop_hbyH0_SI</span> <span class="o">/</span> <span class="nv">coop_Gyr_SI</span></div><div class='line' id='LC103'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_chbyH0_Mpc</span> <span class="o">=</span> <span class="nv">coop_chbyMpcH0</span>  <span class="c">!!alias</span></div><div class='line' id='LC104'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_chbyH0_SI</span> <span class="o">=</span> <span class="nv">coop_chbyH0_Mpc</span> <span class="o">*</span> <span class="nv">coop_Mpc_SI</span>  <span class="c">!!alias</span></div><div class='line' id='LC105'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_H0byh_SI</span> <span class="o">=</span> <span class="mf">1.</span><span class="o">/</span><span class="nv">coop_hbyH0_SI</span></div><div class='line' id='LC106'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_H0byh_eV</span> <span class="o">=</span> <span class="nv">coop_hbar</span><span class="o">/</span><span class="nv">coop_eV_SI</span> <span class="o">*</span> <span class="nv">coop_H0byh_SI</span></div><div class='line' id='LC107'><br/></div><div class='line' id='LC108'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_blackbody_alpha_SI</span> <span class="o">=</span> <span class="p">(</span><span class="nv">coop_pi2</span><span class="o">/</span><span class="mi">3</span><span class="mf">0.</span><span class="err">_</span><span class="nv">dl</span><span class="p">)</span><span class="o">*</span><span class="p">(</span><span class="nv">coop_boltzmann</span><span class="p">)</span> <span class="o">/</span> <span class="p">(</span><span class="nv">coop_PlanckTemperature_SI</span> <span class="o">*</span> <span class="nv">coop_PlanckLength_SI</span><span class="p">)</span><span class="o">**</span><span class="mi">3</span>  <span class="c">!!blackbody radiation density = alpha * (total g) * (Temperature in Kelvin) ^ 4  [g=1 for each spin degree of boson, g=7/8 for each spin degree of fermion; for photon g_total = 2; for neutrinos g_total = 7/8 * number of species] (result is in SI unit J/m^3)</span></div><div class='line' id='LC109'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_Stefan_Boltzmann</span> <span class="o">=</span> <span class="mf">5.670400e-8</span><span class="err">_</span><span class="nv">dl</span> <span class="c">!!W/m^2/K^4</span></div><div class='line' id='LC110'><br/></div><div class='line' id='LC111'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_sigma_thomson_SI</span> <span class="o">=</span> <span class="mf">6.6524616e-29</span><span class="err">_</span><span class="nv">dl</span></div><div class='line' id='LC112'><br/></div><div class='line' id='LC113'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span> <span class="k">parameter</span> <span class="kd">::</span> <span class="nv">coop_arad_SI</span> <span class="o">=</span> <span class="p">(</span><span class="nv">coop_pi</span><span class="o">**</span><span class="mi">5</span> <span class="o">*</span> <span class="mf">8.</span><span class="o">/</span><span class="mi">1</span><span class="mf">5.</span> <span class="p">)</span> <span class="o">*</span> <span class="nv">coop_boltzmann</span> <span class="o">*</span> <span class="p">(</span><span class="nv">coop_boltzmann</span><span class="o">/</span><span class="p">(</span><span class="nv">coop_SpeedOfLight_SI</span> <span class="o">*</span> <span class="nv">coop_Planck</span><span class="p">))</span> <span class="o">**</span><span class="mi">3</span>  </div><div class='line' id='LC114'>&nbsp;&nbsp;&nbsp;&nbsp;<span class="c">!7.565914e-16_dl !radiation constant for u=aT^4</span></div><div class='line' id='LC115'>&nbsp;&nbsp;<span class="c">!! = 2.*coop_blackboday_alpha_SI (since photon has two spin dof)</span></div><div class='line' id='LC116'><br/></div><div class='line' id='LC117'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span> <span class="k">parameter</span> <span class="kd">::</span> <span class="nv">coop_ComptonCT</span> <span class="o">=</span> <span class="p">(</span><span class="mf">8.</span><span class="nv">d0</span><span class="o">/</span><span class="mf">3.</span><span class="nv">d0</span><span class="p">)</span> <span class="o">*</span> <span class="nv">coop_sigma_thomson_SI</span><span class="o">/</span><span class="p">(</span><span class="nv">coop_electron_mass_SI</span> <span class="o">*</span> <span class="nv">coop_SpeedOfLight_SI</span><span class="o">**</span><span class="mi">2</span><span class="p">)</span> <span class="o">*</span> <span class="nv">coop_arad_SI</span>  <span class="o">*</span> <span class="nv">coop_chbyH0_SI</span>  <span class="c">!! ch/H_0 * 8/3 alpha * sigma_T / (m_e c^2)</span></div><div class='line' id='LC118'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span><span class="nv">coop_barssc0</span> <span class="o">=</span> <span class="nv">coop_boltzmann</span> <span class="o">/</span> <span class="nv">coop_hydrogen_mass_SI</span> <span class="o">/</span> <span class="nv">coop_SpeedOfLight_SI</span> <span class="o">**</span> <span class="mi">2</span></div><div class='line' id='LC119'><br/></div><div class='line' id='LC120'><br/></div><div class='line' id='LC121'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span> <span class="nv">coop_rhocritbyh2_MsunbyMpc3</span> <span class="o">=</span> <span class="mf">3.</span><span class="err">_</span><span class="nv">dl</span><span class="o">*</span><span class="p">(</span><span class="nv">coop_PlanckMass_SI</span><span class="o">/</span><span class="nv">coop_Msun_SI</span><span class="p">)</span><span class="o">*</span><span class="p">(</span><span class="nv">coop_Mpc_SI</span><span class="o">/</span><span class="nv">coop_PlanckLength_SI</span><span class="p">)</span> <span class="o">*</span> <span class="p">(</span><span class="mf">1.e5</span><span class="err">_</span><span class="nv">dl</span><span class="o">/</span><span class="nv">coop_SpeedOfLight_SI</span><span class="p">)</span><span class="o">**</span><span class="mi">2</span> <span class="o">/</span> <span class="nv">coop_8pi</span> <span class="c">!!2.77467e11</span></div><div class='line' id='LC122'>&nbsp;&nbsp;<span class="kt">real</span><span class="p">(</span><span class="nv">dl</span><span class="p">),</span><span class="k">parameter</span><span class="kd">::</span><span class="nv">coop_deltac_EdS</span> <span class="o">=</span> <span class="p">(</span><span class="mf">3.</span><span class="err">_</span><span class="nv">dl</span><span class="o">/</span><span class="mf">5.</span><span class="err">_</span><span class="nv">dl</span><span class="p">)</span><span class="o">*</span><span class="p">(</span><span class="mf">3.</span><span class="err">_</span><span class="nv">dl</span><span class="o">/</span><span class="mf">2.</span><span class="err">_</span><span class="nv">dl</span><span class="o">*</span><span class="nv">coop_pi</span><span class="p">)</span><span class="o">**</span><span class="p">(</span><span class="mf">2.</span><span class="err">_</span><span class="nv">dl</span><span class="o">/</span><span class="mf">3.</span><span class="err">_</span><span class="nv">dl</span><span class="p">)</span> <span class="c">!!spherical collapse critical density contrast</span></div><div class='line' id='LC123'><br/></div><div class='line' id='LC124'><br/></div><div class='line' id='LC125'>&nbsp;&nbsp;</div><div class='line' id='LC126'><br/></div><div class='line' id='LC127'><span class="k">end module </span><span class="nv">coop_constants</span></div></pre></div></td>
          </tr>
        </table>
  </div>

  </div>
</div>

<a href="#jump-to-line" rel="facebox[.linejump]" data-hotkey="l" class="js-jump-to-line" style="display:none">Jump to Line</a>
<div id="jump-to-line" style="display:none">
  <form accept-charset="UTF-8" class="js-jump-to-line-form">
    <input class="linejump-input js-jump-to-line-field" type="text" placeholder="Jump to line&hellip;" autofocus>
    <button type="submit" class="button">Go</button>
  </form>
</div>

        </div>

      </div><!-- /.repo-container -->
      <div class="modal-backdrop"></div>
    </div><!-- /.container -->
  </div><!-- /.site -->


    </div><!-- /.wrapper -->

      <div class="container">
  <div class="site-footer">
    <ul class="site-footer-links right">
      <li><a href="https://status.github.com/">Status</a></li>
      <li><a href="http://developer.github.com">API</a></li>
      <li><a href="http://training.github.com">Training</a></li>
      <li><a href="http://shop.github.com">Shop</a></li>
      <li><a href="/blog">Blog</a></li>
      <li><a href="/about">About</a></li>

    </ul>

    <a href="/">
      <span class="mega-octicon octicon-mark-github" title="GitHub"></span>
    </a>

    <ul class="site-footer-links">
      <li>&copy; 2014 <span title="0.05598s from github-fe127-cp1-prd.iad.github.net">GitHub</span>, Inc.</li>
        <li><a href="/site/terms">Terms</a></li>
        <li><a href="/site/privacy">Privacy</a></li>
        <li><a href="/security">Security</a></li>
        <li><a href="/contact">Contact</a></li>
    </ul>
  </div><!-- /.site-footer -->
</div><!-- /.container -->


    <div class="fullscreen-overlay js-fullscreen-overlay" id="fullscreen_overlay">
  <div class="fullscreen-container js-fullscreen-container">
    <div class="textarea-wrap">
      <textarea name="fullscreen-contents" id="fullscreen-contents" class="fullscreen-contents js-fullscreen-contents" placeholder="" data-suggester="fullscreen_suggester"></textarea>
    </div>
  </div>
  <div class="fullscreen-sidebar">
    <a href="#" class="exit-fullscreen js-exit-fullscreen tooltipped tooltipped-w" aria-label="Exit Zen Mode">
      <span class="mega-octicon octicon-screen-normal"></span>
    </a>
    <a href="#" class="theme-switcher js-theme-switcher tooltipped tooltipped-w"
      aria-label="Switch themes">
      <span class="octicon octicon-color-mode"></span>
    </a>
  </div>
</div>



    <div id="ajax-error-message" class="flash flash-error">
      <span class="octicon octicon-alert"></span>
      <a href="#" class="octicon octicon-x close js-ajax-error-dismiss"></a>
      Something went wrong with that request. Please try again.
    </div>


      <script crossorigin="anonymous" src="https://assets-cdn.github.com/assets/frameworks-5bef6dacd990ce272ec009917ceea0b9d96f84b7.js" type="text/javascript"></script>
      <script async="async" crossorigin="anonymous" src="https://assets-cdn.github.com/assets/github-1b13ad871c0387dc91be6a720577ab10d28b22af.js" type="text/javascript"></script>
      
      
  </body>
</html>

