package main

import (
	"fmt"
	"path"
	"path/filepath"
	"strings"

	"github.com/jjtimmons/repp/cmd"
	"github.com/spf13/cobra/doc"
)

// https://pmarsceill.github.io/just-the-docs/docs/navigation-structure/
const rootCmd = `---
layout: default
title: %s
nav_order: %d
has_children: true
permalink: /
---
`

// child command without children
const childCmd = `---
layout: default
title: %s
parent: %s
nav_order: %d
---
`

// child with children
const childParentCmd = `---
layout: default
title: %s
parent: %s
nav_order: %d
has_children: true
---
`

// grandchildren
const grandchildCmd = `---
layout: default
title: %s
parent: %s
grand_parent: %s
nav_order: %d
---
`

// docType codes whether the command is a grandchild, child, etc
type docType int

const (
	root docType = iota
	child
	childParent
	grandchild
)

// meta is for describing the position/info for a command doc page
type meta struct {
	docType     docType
	title       string
	navOrder    int
	hasChildren bool
	parent      string
	grandParent string
}

// map from the base Markdown file name to its build meta
var metaMap = map[string]meta{
	"repp": meta{
		root,
		"repp",
		0,
		true,
		"",
		"",
	},
	"repp_make": meta{
		childParent,
		"make",
		0,
		true,
		"repp",
		"",
	},
	"repp_make_sequence": meta{
		grandchild,
		"sequence",
		0,
		false,
		"make",
		"repp",
	},
	"repp_make_features": meta{
		grandchild,
		"features",
		1,
		false,
		"make",
		"repp",
	},
	"repp_make_fragments": meta{
		grandchild,
		"fragments",
		2,
		false,
		"make",
		"repp",
	},
	"repp_find": meta{
		childParent,
		"find",
		1,
		true,
		"repp",
		"",
	},
	"repp_find_sequence": meta{
		grandchild,
		"sequence",
		0,
		false,
		"find",
		"repp",
	},
	"repp_find_fragment": meta{
		grandchild,
		"fragment",
		1,
		false,
		"find",
		"repp",
	},
	"repp_find_feature": meta{
		grandchild,
		"feature",
		2,
		false,
		"find",
		"repp",
	},
	"repp_find_enzyme": meta{
		grandchild,
		"enzyme",
		3,
		false,
		"find",
		"repp",
	},
	"repp_set": meta{
		childParent,
		"set",
		2,
		true,
		"repp",
		"",
	},
	"repp_set_feature": meta{
		grandchild,
		"feature",
		0,
		false,
		"set",
		"repp",
	},
	"repp_set_enzyme": meta{
		grandchild,
		"enzyme",
		1,
		false,
		"set",
		"repp",
	},
	"repp_delete": meta{
		childParent,
		"delete",
		3,
		true,
		"repp",
		"",
	},
	"repp_delete_feature": meta{
		grandchild,
		"feature",
		0,
		false,
		"delete",
		"repp",
	},
	"repp_annotate": meta{
		child,
		"annotate",
		4,
		false,
		"repp",
		"",
	},
}

// makeDocs parses the custom commands and outputs Markdown documentation files
func makeDocs() {
	if err := doc.GenMarkdownTreeCustom(cmd.RootCmd, "./docs", filePrepender, linkHandler); err != nil {
		fmt.Println(err.Error())
	}
}

// filePrepender adds YAML headings that are required by the just-the-docs theme
// https://github.com/spf13/cobra/blob/master/doc/md_docs.md
// https://pmarsceill.github.io/just-the-docs/docs/navigation-structure/
func filePrepender(filename string) string {
	name := filepath.Base(filename)
	base := strings.TrimSuffix(name, path.Ext(name))
	m, _ := metaMap[base]

	switch m.docType {
	case root:
		return fmt.Sprintf(rootCmd, m.title, m.navOrder)
	case child:
		return fmt.Sprintf(childCmd, m.title, m.parent, m.navOrder)
	case childParent:
		return fmt.Sprintf(childParentCmd, m.title, m.parent, m.navOrder)
	case grandchild:
		return fmt.Sprintf(grandchildCmd, m.title, m.parent, m.grandParent, m.navOrder)
	}

	return ""
}

/// linkHandler returns the URL to a documentation page
func linkHandler(filename string) string {
	name := filepath.Base(filename)
	base := strings.TrimSuffix(name, path.Ext(name))

	if base == "repp" {
		return "/"
	}
	return base
}
