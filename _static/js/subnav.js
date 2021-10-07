// Sort toctree items by class
// in order to correctly display nav
const subNav = document.getElementById('subnav')
sortList(subNav.children[0].children[0])
subNav.children[0].style.display = 'block'

function sortList(list, level) {
   if (level === undefined)
      level = 0

   if (!list.children)
      return

   const length = list.children.length
   if (length === 0)
      return

   for (let i = 0; i < length; i++)
      sortItem(list.children[i], level)
}

function sortItem(item, level) {
   if (!item.children)
      return

   const length = item.children.length
   if (length === 0)
      return

   const link = item.children[0]
   let ref = item.children[0].href
   ref = ref.replace(/^.*[\\\/]/, '')

   if (ref.search(/#./) != -1) {
      link.classList.add('anchor')
      link.parentNode.classList.add('anchor')
      // link.parentNode.parentNode.classList.add('anchor')
   }
   else {
      link.classList.add('link', 'level'+level)
      link.parentNode.classList.add('link', 'level'+level)
      // link.parentNode.parentNode.classList.add('link', 'level'+level)
   }

   level++
   for (let i = 1; i < length; i++)
      sortList(item.children[i], level)
}
