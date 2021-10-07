function openNav() {
   document.getElementById("sideNav").style.width = "300px";
}

function closeNav() {
   document.getElementById("sideNav").style.width = "0";
}

function previousNav(item, level) {
   showSideNav(item.parent, level-1)
}

function nextNav(e, item, level) {
   e.preventDefault()
   e.stopPropagation()
   showSideNav(item, level+1)
}

// Convert list into div
const toctree = nav.children[0].children[0]

const tree = []
let currentLevel = 0
let currentItem = undefined
let currentParent = undefined
toctree2tree(toctree, 0, undefined)

function toctree2tree(toctree, level, parent) {
   //toctree === ul
   const length = toctree.children.length

   for (let i = 0; i < length; i++) {
      const item = toctree.children[i] //item === li
      const link = item.children[0] //link === a

      const object = {
         'level': level,
         'parent': parent,
         'className': 'level'+level+' nav-item',
         'innerHTML': link.innerHTML,
         'href': link.href,
         'children': []
      }

      if (item.classList.contains('current')) {
         object.className += ' current'
         if (level > currentLevel) {
            currentItem = object
            currentLevel = level
            currentParent = parent
         }
      }

      let ref = link.href
      ref = ref.replace(/^.*[\\\/]/, '')

      if (ref.search(/#./) != -1)
         object.type = 'anchor'
      else
         object.type = 'link'

      tree.push(object)

      if (item.children[1])   //=== ul
         toctree2tree(item.children[1], level+1, object)

      if (parent)
         parent.children.push(object)
   }
}

// Remove toctree
nav.removeChild(nav.children[0])

// Append level0 to staticNav
const staticNav = document.getElementById('staticNav')
tree.forEach(function(item) {
   if (item.level === 0 && item.innerHTML !== 'References') {
      const div = document.createElement('div')
      div.className = item.className
      div.onclick = function() { window.location.href=item.href }
      div.innerHTML = item.innerHTML
      staticNav.appendChild(div)
   }
})

// Append current level to sideNav
if (!currentItem)
   showSideNav(currentParent, currentLevel)
else if (currentItem && currentItem.children && currentItem.children.length === 0)
   showSideNav(currentParent, currentLevel)
else
   showSideNav(currentItem, currentLevel+1)

function showSideNav(parent, level) {
   const dsideNav = document.getElementById('dynamicSideNav')

   dsideNav.innerHTML = ''

   const sideNavPrevious = document.createElement('a')
   sideNavPrevious.href = "#"
   sideNavPrevious.className = 'previousNav'
   sideNavPrevious.innerHTML = '<i class="fas fa-arrow-left"></i>'

   if (parent) { // Show current title & show previous nav button
      sideNavPrevious.onclick = function() { previousNav(parent, level) }
      dsideNav.appendChild(sideNavPrevious)

      const sideNavTitle = document.createElement('div')
      sideNavTitle.className = 'sidenav-title'
      sideNavTitle.innerHTML = parent.innerHTML
      dsideNav.appendChild(sideNavTitle)
   }

   const globalDiv = document.createElement('div')
   globalDiv.className = 'nav-items'

   tree.forEach(function(item) {
   if ((item.level === level) && (item.parent === parent)) {
         const div = document.createElement('div')
         div.className = item.className
         div.onclick = function() { window.location.href=item.href; closeNav() }

         const p = document.createElement('p')
         p.innerHTML = item.innerHTML

         div.appendChild(p)

         if (item.children.length > 0) {
            const next = document.createElement('a')
            next.className = 'sidenav-next'
            next.href = "#"
            next.onclick = function(e) { nextNav(e, item, level) }

            if (item.children[0].type === 'link')
               next.innerHTML = '<i class="fas fa-chevron-right"></i>'
            else
               next.innerHTML = '<i class="fas fa-list-ul"></i>'

            div.appendChild(next)
         }

         globalDiv.appendChild(div)
      }
   })
   dsideNav.appendChild(globalDiv)
}
